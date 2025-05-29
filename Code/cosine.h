#ifndef COSINE_H
#define COSINE_H

#include <math.h>
#include <stdint.h>

// Simple struct for vector representation
// K-mer is now a 64-bit unsigned integer
typedef struct {
  uint64_t kmer; // Changed from char*
  int count;
} KmerCount;

// Forward declaration of KmerTable from kmer_counter.h
typedef struct KmerTable {
  size_t size;     // Number of unique k-mers currently in the table
  size_t capacity; // Total slots in the hash table
  KmerCount *data; // Array of KmerCount structs (hash table slots)
} KmerTable;

int get_count_from_table(KmerTable *table, uint64_t kmer);
double cosine_similarity(KmerTable *v1_table, KmerTable *v2_table);

inline double cosine_similarity(KmerTable *v1_table, KmerTable *v2_table) {
  double dot = 0.0;
  double norm1 = 0.0;
  double norm2 = 0.0;

  // Iterate through the first table's unique k-mers
  // This assumes KmerTable has an accessible data array for iteration.
  // The KmerTable's 'data' array contains the (kmer, count) pairs.
  for (size_t i = 0; i < v1_table->capacity; i++) {
    if (v1_table->data[i].kmer != 0) { // Check if slot is occupied (assuming 0
                                       // is never a valid k-mer hash)
      uint64_t current_kmer_v1 = v1_table->data[i].kmer;
      int count1 = v1_table->data[i].count;
      int count2 =
          get_count_from_table(v2_table, current_kmer_v1); // Hash lookup in v2

      dot += (double)count1 * count2;
      norm1 += (double)count1 * count1;
    }
  }

  // Now calculate norm2 by iterating through the second table
  for (size_t i = 0; i < v2_table->capacity; i++) {
    if (v2_table->data[i].kmer != 0) { // Check if slot is occupied
      int count2 = v2_table->data[i].count;
      norm2 += (double)count2 * count2;
    }
  }

  if (norm1 == 0.0 || norm2 == 0.0)
    return 0.0;
  return dot / (sqrt(norm1) * sqrt(norm2));
}

#endif // COSINE_H