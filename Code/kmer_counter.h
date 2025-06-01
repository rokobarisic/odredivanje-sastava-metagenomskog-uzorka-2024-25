#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <ctype.h>
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_KMER_LEN 31
#define INITIAL_TABLE_CAPACITY 8192
#define EMPTY_KMER UINT64_MAX

// Pre-computed lookup table for base conversion (much faster than switch)
// Compiler hack, set all 256 values to 4, but then immediately set some of them
// to something else. It looks redundant but works like a charm.
static const uint64_t BASE_LOOKUP[256] = {
    [0 ... 255] = 4, // Invalid base marker
    ['A'] = 0,       ['a'] = 0, ['C'] = 1, ['c'] = 1,
    ['G'] = 2,       ['g'] = 2, ['T'] = 3, ['t'] = 3};

typedef struct {
  uint64_t kmer;
  uint32_t count;
  uint16_t distance;
} RobinHoodEntry;

typedef struct {
  size_t size;
  size_t capacity;
  RobinHoodEntry *data;
  uint32_t k;
  uint64_t mask;
} RobinHoodTable;

typedef struct {
  RobinHoodTable **tables;
  double *similarities;
  size_t count;
} BatchResult;

static inline uint64_t base_to_bits_fast(char c) {
  return BASE_LOOKUP[(unsigned char)c];
}

static inline uint64_t hash_kmer_fast(uint64_t kmer) {
  kmer ^= kmer >> 32;
  kmer *= 0x9e3779b97f4a7c15ULL;
  kmer ^= kmer >> 25;
  kmer *= 0x9e3779b97f4a7c15ULL;
  kmer ^= kmer >> 16;
  return kmer;
}

static RobinHoodTable *create_robin_hood_table(int k) {
  RobinHoodTable *table = (RobinHoodTable *)malloc(sizeof(RobinHoodTable));
  if (!table)
    return NULL;

  table->size = 0;
  table->capacity = INITIAL_TABLE_CAPACITY;
  table->k = k;
  table->mask = (1ULL << (2 * k)) - 1;
  table->data = (RobinHoodEntry *)aligned_alloc(64, table->capacity *
                                                        sizeof(RobinHoodEntry));
  if (!table->data) {
    free(table);
    return NULL;
  }

  for (size_t i = 0; i < table->capacity; ++i) {
    table->data[i].kmer = EMPTY_KMER;
  }

  return table;
}

static void reinsert_robin_hood(RobinHoodTable *table, uint64_t kmer,
                                uint32_t count) {
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
      RobinHoodEntry temp = *slot;
      entry.distance = distance;
      *slot = entry;
      entry = temp;
      distance = entry.distance;
    }
    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

static void resize_robin_hood_table(RobinHoodTable *table) {
  size_t old_capacity = table->capacity;
  RobinHoodEntry *old_data = table->data;

  table->capacity *= 2;
  table->size = 0;
  table->data = (RobinHoodEntry *)aligned_alloc(64, table->capacity *
                                                        sizeof(RobinHoodEntry));
  if (!table->data) {
    table->capacity = old_capacity;
    table->data = old_data;
    return;
  }

  for (size_t i = 0; i < table->capacity; ++i)
    table->data[i].kmer = EMPTY_KMER;

  for (size_t i = 0; i < old_capacity; i++) {
    if (old_data[i].kmer != EMPTY_KMER) {
      reinsert_robin_hood(table, old_data[i].kmer, old_data[i].count);
    }
  }

  free(old_data);
}

static void insert_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  if (table->size * 10 > table->capacity * 7) {
    resize_robin_hood_table(table);
  }

  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1);
  uint16_t distance = 0;
  RobinHoodEntry entry = {kmer, 1, 0};

  while (1) {
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
      RobinHoodEntry temp = *slot;
      entry.distance = distance;
      *slot = entry;
      entry = temp;
      distance = entry.distance;
    }

    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

static RobinHoodTable *count_kmers_optimized(const char *sequence, int k) {
  if (k <= 0 || k > MAX_KMER_LEN)
    return NULL;

  RobinHoodTable *table = create_robin_hood_table(k);
  if (!table)
    return NULL;

  size_t len = strlen(sequence);
  if (len < k)
    return table;

  uint64_t current_kmer = 0;
  uint64_t mask = table->mask;
  int valid_bases = 0;

  for (size_t i = 0; i < len; i++) {
    uint64_t base_bits = base_to_bits_fast(sequence[i]);

    if (base_bits == 4) {
      current_kmer = 0;
      valid_bases = 0;
      continue;
    }

    current_kmer = ((current_kmer << 2) | base_bits) & mask;
    valid_bases++;

    if (valid_bases >= k)
      insert_robin_hood(table, current_kmer);
  }

  return table;
}

static uint32_t get_count_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  if (!table || table->size == 0)
    return 0;

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

static double cosine_similarity_optimized(RobinHoodTable *table1,
                                          RobinHoodTable *table2) {
  if (!table1 || !table2 || table1->size == 0 || table2->size == 0)
    return 0.0;

  RobinHoodTable *smaller = (table1->size <= table2->size) ? table1 : table2;
  RobinHoodTable *larger = (table1->size <= table2->size) ? table2 : table1;

  double dot_product = 0.0, norm1_sq = 0.0, norm2_sq = 0.0;

  for (size_t i = 0; i < smaller->capacity; i++) {
    if (smaller->data[i].kmer != EMPTY_KMER) {
      uint32_t count1 = smaller->data[i].count;
      uint32_t count2 = get_count_robin_hood(larger, smaller->data[i].kmer);

      dot_product += (double)count1 * count2;
    }
  }

  for (size_t i = 0; i < table1->capacity; i++) {
    if (table1->data[i].kmer != EMPTY_KMER) {
      double count = (double)table1->data[i].count;
      norm1_sq += count * count;
    }
  }

  for (size_t i = 0; i < table2->capacity; i++) {
    if (table2->data[i].kmer != EMPTY_KMER) {
      double count = (double)table2->data[i].count;
      norm2_sq += count * count;
    }
  }

  if (norm1_sq == 0.0 || norm2_sq == 0.0)
    return 0.0;

  return dot_product / (sqrt(norm1_sq) * sqrt(norm2_sq));
}

static void free_robin_hood_table(RobinHoodTable *table) {
  if (!table)
    return;
  free(table->data);
  free(table);
}

#endif // KMER_COUNTER_H
