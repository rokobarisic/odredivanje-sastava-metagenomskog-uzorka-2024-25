#include <x86intrin.h> // For _mm_prefetch()

#include "kmer_counter.h"

// Pre-computed lookup table for faster base conversion
const uint64_t BASE_LOOKUP[256] = {
    [0 ... 255] = 4, ['A'] = 0, ['a'] = 0, ['C'] = 1, ['c'] = 1,
    ['G'] = 2,       ['g'] = 2, ['T'] = 3, ['t'] = 3};

inline uint64_t base_to_bits_fast(char c) {
  return BASE_LOOKUP[(unsigned char)c];
}

inline uint64_t hash_kmer_fast(uint64_t kmer) {
  return (kmer * 0x9e3779b97f4a7c13ULL) >> 16;
}

inline double fast_inv_sqrt(double x) {
  union {
    double d;
    uint64_t i;
  } conv = {x};
  conv.i = 0x5fe6ec85e7de30daULL - (conv.i >> 1);
  double y = conv.d;
  return y * (1.5 - 0.5 * x * y * y); // 1 Newton-Raphson iteration
}

RobinHoodTable *create_robin_hood_table(int k) {
  RobinHoodTable *table = aligned_alloc(64, sizeof(RobinHoodTable));
  if (!table)
    return NULL;

  table->size = 0;
  table->capacity = INITIAL_TABLE_CAPACITY;
  table->k = k;
  table->mask = (1ULL << (2 * k)) - 1;
  table->data = aligned_alloc(64, table->capacity * sizeof(RobinHoodEntry));
  if (!table->data) {
    free(table);
    return NULL;
  }

  for (size_t i = 0; i < table->capacity; ++i)
    table->data[i].kmer = EMPTY_KMER;

  return table;
}

void insert_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  if ((table->size << 1) > table->capacity)
    resize_robin_hood_table(table);

  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1);
  uint16_t distance = 0;
  RobinHoodEntry entry = {kmer, 1, 0};

  while (1) {
    _mm_prefetch(&table->data[(pos + 1) & (table->capacity - 1)], _MM_HINT_T0);
    RobinHoodEntry *slot = &table->data[pos];

    if (slot->kmer == EMPTY_KMER) {
      entry.distance = distance;
      *slot = entry;
      table->size++;
      return;
    }

    if (slot->kmer == kmer) {
      slot->count++;
      return;
    }

    if (distance > slot->distance) {
      RobinHoodEntry tmp = *slot;
      entry.distance = distance;
      *slot = entry;
      entry = tmp;
      distance = entry.distance;
    }

    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

RobinHoodTable *cnt_kmer(const char *sequence, int k) {
  if (k <= 0 || k > MAX_KMER_LEN)
    return NULL;

  RobinHoodTable *table = create_robin_hood_table(k);
  if (!table)
    return NULL;

  size_t len = strlen(sequence);
  if (len < (size_t)k)
    return table;

  uint64_t current_kmer = 0;
  uint64_t mask = table->mask;
  int valid_bases = 0;

  for (size_t i = 0; i < len; i++) {
    uint64_t base = base_to_bits_fast(sequence[i]);
    if (base > 3) {
      valid_bases = 0;
      current_kmer = 0;
      continue;
    }

    current_kmer = ((current_kmer << 2) | base) & mask;
    valid_bases++;

    if (valid_bases >= k)
      insert_robin_hood(table, current_kmer);

    _mm_prefetch(&sequence[i + 16], _MM_HINT_T0);
  }

  return table;
}

uint32_t get_count_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1);
  uint16_t distance = 0;

  while (1) {
    RobinHoodEntry *slot = &table->data[pos];
    if (slot->kmer == EMPTY_KMER || distance > slot->distance)
      return 0;
    if (slot->kmer == kmer)
      return slot->count;

    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

double cos_similarity(RobinHoodTable *table1, RobinHoodTable *table2) {
  if (!table1 || !table2 || table1->size == 0 || table2->size == 0)
    return 0.0;

  double dot_product = 0.0;
  double norm1_sq = 0.0;
  double norm2_sq = 0.0;

  for (size_t i = 0; i < table1->capacity; i++) {
    RobinHoodEntry *e = &table1->data[i];
    if (e->kmer != EMPTY_KMER) {
      double count1 = (double)e->count;
      norm1_sq += count1 * count1;
      dot_product += count1 * get_count_robin_hood(table2, e->kmer);
    }
  }

  for (size_t i = 0; i < table2->capacity; i++) {
    RobinHoodEntry *e = &table2->data[i];
    if (e->kmer != EMPTY_KMER) {
      double count2 = (double)e->count;
      norm2_sq += count2 * count2;
    }
  }

  if (norm1_sq == 0.0 || norm2_sq == 0.0)
    return 0.0;

  // Use fast inverse sqrt for approximate cosine similarity
  double denom = fast_inv_sqrt(norm1_sq) * fast_inv_sqrt(norm2_sq);
  return dot_product * denom;
}

void resize_robin_hood_table(RobinHoodTable *table) {
  size_t old_capacity = table->capacity;
  RobinHoodEntry *old_data = table->data;

  table->capacity *= 2;
  table->size = 0;
  table->data = aligned_alloc(64, table->capacity * sizeof(RobinHoodEntry));
  if (!table->data) {
    table->capacity = old_capacity;
    table->data = old_data;
    return;
  }

  for (size_t i = 0; i < table->capacity; ++i)
    table->data[i].kmer = EMPTY_KMER;

  for (size_t i = 0; i < old_capacity; ++i) {
    if (old_data[i].kmer != EMPTY_KMER) {
      reinsert_robin_hood(table, old_data[i].kmer, old_data[i].count);
    }
  }

  free(old_data);
}

void reinsert_robin_hood(RobinHoodTable *table, uint64_t kmer, uint32_t count) {
  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1);
  uint16_t distance = 0;
  RobinHoodEntry entry = {kmer, count, distance};

  while (1) {
    RobinHoodEntry *slot = &table->data[pos];
    if (slot->kmer == EMPTY_KMER) {
      *slot = entry;
      table->size++;
      return;
    }

    if (distance > slot->distance) {
      RobinHoodEntry tmp = *slot;
      entry.distance = distance;
      *slot = entry;
      entry = tmp;
      distance = entry.distance;
    }

    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

void free_robin_hood_table(RobinHoodTable *table) {
  if (table) {
    free(table->data);
    free(table);
  }
}
