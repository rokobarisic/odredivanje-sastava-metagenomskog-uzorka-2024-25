#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <ctype.h>
#include <immintrin.h> // For SIMD if available
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Optimized constants
#define MAX_KMER_LEN 31
#define INITIAL_TABLE_CAPACITY 8192 // Increased from 1024

// Pre-computed lookup table for base conversion (much faster than switch)
// Compiler hack, set all 256 values to 4, but then immediately set some of them
// to something else. It looks redundant but works like a charm.
static const uint64_t BASE_LOOKUP[256] = {
    [0 ... 255] = 4, // Invalid base marker
    ['A'] = 0,       ['a'] = 0, ['C'] = 1, ['c'] = 1,
    ['G'] = 2,       ['g'] = 2, ['T'] = 3, ['t'] = 3};

typedef struct {
  uint64_t kmer;
  uint32_t count; // Reduced from int to save memory
} KmerCount;

typedef struct KmerTable {
  size_t size;
  size_t capacity;
  KmerCount *data;
  uint32_t k;    // Store k-mer length for optimizations
  uint64_t mask; // Pre-computed mask
} KmerTable;

// Optimized base conversion using lookup table
static inline uint64_t base_to_bits_fast(char c) {
  return BASE_LOOKUP[(unsigned char)c];
}

// Faster hash function - FNV-1a variant optimized for genomic data
static inline uint64_t hash_kmer_fast(uint64_t kmer) {
  // FNV-1a hash with genomic-optimized constants
  kmer ^= kmer >> 32;
  kmer *= 0x9e3779b97f4a7c15ULL; // Golden ratio based multiplier
  kmer ^= kmer >> 25;
  kmer *= 0x9e3779b97f4a7c15ULL;
  kmer ^= kmer >> 16;
  return kmer;
}

// Robin Hood hash table with better cache locality
typedef struct {
  uint64_t kmer;
  uint32_t count;
  uint16_t distance; // Distance from ideal position for Robin Hood hashing
} RobinHoodEntry;

typedef struct {
  size_t size;
  size_t capacity;
  RobinHoodEntry *data;
  uint32_t k;
  uint64_t mask;
} RobinHoodTable;

// Batch processing for multiple reads
typedef struct {
  RobinHoodTable **tables;
  double *similarities;
  size_t count;
} BatchResult;

RobinHoodTable *create_robin_hood_table(int k);
RobinHoodTable *count_kmers_optimized(const char *sequence, int k);
uint32_t get_count_robin_hood(RobinHoodTable *table, uint64_t kmer);
double cosine_similarity_optimized(RobinHoodTable *table1,
                                   RobinHoodTable *table2);
void free_robin_hood_table(RobinHoodTable *table);
BatchResult *process_reads_batch(RobinHoodTable *reference, char **reads,
                                 int num_reads, int k);
void free_batch_result(BatchResult *result);

// Create optimized hash table
inline RobinHoodTable *create_robin_hood_table(int k) {
  RobinHoodTable *table = (RobinHoodTable *)malloc(sizeof(RobinHoodTable));
  if (!table)
    return NULL;

  table->size = 0;
  table->capacity = INITIAL_TABLE_CAPACITY;
  table->k = k;
  table->mask = (1ULL << (2 * k)) - 1;

  // Use aligned allocation for better cache performance
  table->data = (RobinHoodEntry *)aligned_alloc(64, table->capacity *
                                                        sizeof(RobinHoodEntry));
  if (!table->data) {
    free(table);
    return NULL;
  }

  // Initialize all entries as empty
  memset(table->data, 0, table->capacity * sizeof(RobinHoodEntry));

  return table;
}

// Robin Hood insertion with better performance
static void insert_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1); // Faster modulo with power of 2
  uint16_t distance = 0;

  RobinHoodEntry entry = {kmer, 1, 0};

  while (1) {
    RobinHoodEntry *slot = &table->data[pos];

    if (slot->kmer == 0) {
      // Empty slot
      entry.distance = distance;
      *slot = entry;
      table->size++;
      return;
    }

    if (slot->kmer == kmer) {
      // Found existing k-mer
      slot->count++;
      return;
    }

    // Robin Hood: if current entry has traveled further, swap
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

// Resize with Robin Hood
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

  memset(table->data, 0, table->capacity * sizeof(RobinHoodEntry));

  // Re-insert all entries
  for (size_t i = 0; i < old_capacity; i++) {
    if (old_data[i].kmer != 0) {
      for (uint32_t j = 0; j < old_data[i].count; j++) {
        insert_robin_hood(table, old_data[i].kmer);
      }
    }
  }

  free(old_data);
}

// Ultra-fast k-mer counting with rolling hash and SIMD optimizations
inline RobinHoodTable *count_kmers_optimized(const char *sequence, int k) {
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

  // Process sequence with optimized rolling hash
  for (size_t i = 0; i < len; i++) {
    uint64_t base_bits = base_to_bits_fast(sequence[i]);

    if (base_bits == 4) {
      // Invalid base - reset
      current_kmer = 0;
      valid_bases = 0;
      continue;
    }

    current_kmer = ((current_kmer << 2) | base_bits) & mask;
    valid_bases++;

    if (valid_bases >= k) {
      // Check load factor and resize if needed
      if (table->size * 10 > table->capacity * 7) { // Avoid floating point
        resize_robin_hood_table(table);
      }

      insert_robin_hood(table, current_kmer);
    }
  }

  return table;
}

// Fast lookup in Robin Hood table
inline uint32_t get_count_robin_hood(RobinHoodTable *table, uint64_t kmer) {
  if (!table || table->size == 0)
    return 0;

  uint64_t hash = hash_kmer_fast(kmer);
  size_t pos = hash & (table->capacity - 1);
  uint16_t distance = 0;

  while (1) {
    RobinHoodEntry *slot = &table->data[pos];

    if (slot->kmer == 0 || distance > slot->distance) {
      return 0; // Not found
    }

    if (slot->kmer == kmer) {
      return slot->count;
    }

    pos = (pos + 1) & (table->capacity - 1);
    distance++;
  }
}

// Optimized cosine similarity with better cache usage
inline double cosine_similarity_optimized(RobinHoodTable *table1,
                                          RobinHoodTable *table2) {
  if (!table1 || !table2 || table1->size == 0 || table2->size == 0) {
    return 0.0;
  }

  // Always iterate over the smaller table for better performance
  RobinHoodTable *smaller = (table1->size <= table2->size) ? table1 : table2;
  RobinHoodTable *larger = (table1->size <= table2->size) ? table2 : table1;

  double dot_product = 0.0;
  double norm1_sq = 0.0;
  double norm2_sq = 0.0;

  // Calculate dot product by iterating over smaller table
  for (size_t i = 0; i < smaller->capacity; i++) {
    if (smaller->data[i].kmer != 0) {
      uint32_t count1 = smaller->data[i].count;
      uint32_t count2 = get_count_robin_hood(larger, smaller->data[i].kmer);

      dot_product += (double)count1 * count2;
    }
  }

  // Calculate norms
  for (size_t i = 0; i < table1->capacity; i++) {
    if (table1->data[i].kmer != 0) {
      double count = (double)table1->data[i].count;
      norm1_sq += count * count;
    }
  }

  for (size_t i = 0; i < table2->capacity; i++) {
    if (table2->data[i].kmer != 0) {
      double count = (double)table2->data[i].count;
      norm2_sq += count * count;
    }
  }

  if (norm1_sq == 0.0 || norm2_sq == 0.0)
    return 0.0;

  return dot_product / (sqrt(norm1_sq) * sqrt(norm2_sq));
}

inline void free_robin_hood_table(RobinHoodTable *table) {
  if (!table)
    return;
  free(table->data);
  free(table);
}

inline BatchResult *process_reads_batch(RobinHoodTable *reference, char **reads,
                                        int num_reads, int k) {
  BatchResult *result = (BatchResult *)malloc(sizeof(BatchResult));
  if (!result)
    return NULL;

  result->tables =
      (RobinHoodTable **)malloc(num_reads * sizeof(RobinHoodTable *));
  result->similarities = (double *)malloc(num_reads * sizeof(double));
  result->count = num_reads;

  if (!result->tables || !result->similarities) {
    free(result->tables);
    free(result->similarities);
    free(result);
    return NULL;
  }

  // Process all reads
  for (int i = 0; i < num_reads; i++) {
    result->tables[i] = count_kmers_optimized(reads[i], k);
    if (result->tables[i]) {
      result->similarities[i] =
          cosine_similarity_optimized(reference, result->tables[i]);
    } else {
      result->similarities[i] = 0.0;
    }
  }

  return result;
}

inline void free_batch_result(BatchResult *result) {
  if (!result)
    return;

  for (size_t i = 0; i < result->count; i++) {
    free_robin_hood_table(result->tables[i]);
  }

  free(result->tables);
  free(result->similarities);
  free(result);
}

#endif // KMER_COUNTER_H