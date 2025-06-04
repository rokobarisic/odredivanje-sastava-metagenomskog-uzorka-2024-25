#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <ctype.h>
#include <immintrin.h> // For AVX/SSE intrinsics
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @def MAX_KMER_LEN
 * @brief Maximum supported k-mer length.
 *
 * Defines the maximum 'k' value for k-mers that can be processed.
 * Limited by the `uint64_t` representation (32 bases, 2 bits per base).
 */
#define MAX_KMER_LEN 31

/**
 * @def INITIAL_TABLE_CAPACITY
 * @brief Initial capacity for the Robin Hood hash table.
 *
 * Initial number of slots allocated for the hash table.
 * This value should ideally be a power of 2 for efficient modulo operations.
 */
#define INITIAL_TABLE_CAPACITY 8192

/**
 * @def EMPTY_KMER
 * @brief Sentinel value to denote an empty slot in the Robin Hood hash table.
 *
 * This value is used to mark an unoccupied or deleted slot in the hash table,
 * leveraging the maximum possible `uint64_t` value which is unlikely to be
 * a valid k-mer representation.
 */
#define EMPTY_KMER UINT64_MAX

/**
 * @struct RobinHoodEntry
 * @brief Represents a single entry within the Robin Hood hash table.
 *
 * Each entry stores a k-mer (encoded as a 64-bit integer), its occurrence
 * count, and its 'distance' from its ideal hash position, which is crucial
 * for the Robin Hood hashing strategy.
 */
typedef struct {
  uint64_t kmer;     /**< @brief The k-mer sequence encoded as a 64-bit unsigned
                        integer. */
  uint32_t count;    /**< @brief The frequency count of this k-mer. */
  uint16_t distance; /**< @brief The "probe distance" or how far this entry is
                        from its preferred slot. */
} RobinHoodEntry;

/**
 * @struct RobinHoodTable
 * @brief Represents the Robin Hood hash table structure.
 *
 * Encapsulates the state of the hash table, including its current
 * number of elements, total capacity, the array of entries,
 * the k-mer length it's configured for and a bitmask used for efficient
 * k-mer hashing.
 */
typedef struct {
  size_t
      size; /**< @brief Current number of active k-mer entries in the table. */
  size_t capacity; /**< @brief Total number of slots (entries) in the table. */
  RobinHoodEntry *data; /**< @brief Pointer to the dynamically allocated array
                           of RobinHoodEntry. */
  uint32_t k; /**< @brief The k-mer length this table is designed to store. */
  uint64_t
      mask; /**< @brief A bitmask used for efficient k-mer encoding/decoding. */
} RobinHoodTable;

/**
 * @brief Converts a DNA base character to its 2-bit integer representation.
 *
 * This function takes a single character representing a DNA base (A, C, G, T/U)
 * and converts it into its corresponding 2-bit integer value. It's
 * case-insensitive.
 *
 * @param c The DNA base character ('A', 'C', 'G', 'T', 'U').
 * @return The 2-bit integer representation (A=00, C=01, G=10, T/U=11).
 */
uint64_t base_to_bits_fast(char c);

/**
 * @brief Computes the hash value for a given k-mer (represented as a 64-bit
 * integer).
 *
 * Calculates the initial hash index for a k-mer within the Robin Hood
 * hash table. The specific hashing algorithm is optimized for speed.
 *
 * @param kmer The k-mer encoded as `uint64_t`.
 * @return The initial hash index (slot) for the k-mer.
 */
uint64_t hash_kmer_fast(uint64_t kmer);

/**
 * @brief Creates and initializes a new Robin Hood hash table.
 *
 * Allocates memory for a new hash table structure and its underlying data
 * array. Initializes all slots as empty and sets up the table's properties
 * based on the provided k-mer length.
 *
 * @param k The k-mer length for which this table will be used.
 * @return A pointer to the newly created and initialized RobinHoodTable, or
 * NULL on memory allocation failure.
 */
RobinHoodTable *create_robin_hood_table(int k);

/**
 * @brief Reinserts k-mer and its count into the Robin Hood hash table during
 * resizing.
 *
 * Internal helper function used specifically when the hash table is
 * resized. Handles the placement of an existing k-mer and its count into the
 * new, larger table following the Robin Hood hashing strategy.
 *
 * @param table Pointer to the RobinHoodTable where the k-mer will be
 * reinserted.
 * @param kmer The k-mer (encoded as `uint64_t`) to reinsert.
 * @param count The frequency count of the k-mer to reinsert.
 */
void reinsert_robin_hood(RobinHoodTable *table, uint64_t kmer, uint32_t count);

/**
 * @brief Resizes the Robin Hood hash table when it becomes too full.
 *
 * Doubles the capacity of the hash table and rehashes all existing entries
 * into the new, larger table. This operation is critical for maintaining
 * performance as the table grows.
 *
 * @param table A pointer to the RobinHoodTable to be resized.
 */
void resize_robin_hood_table(RobinHoodTable *table);

/**
 * @brief Inserts a new k-mer into the Robin Hood hash table or increments its
 * count if it exists.
 *
 * This function handles the core logic of inserting a k-mer. If the k-mer is
 * new, it's added with a count of 1. If it already exists, its count is
 * incremented. It applies the Robin Hood displacement strategy to handle
 * collisions.
 *
 * @param table A pointer to the RobinHoodTable where the k-mer will be
 * inserted.
 * @param kmer The k-mer (encoded as `uint64_t`) to insert.
 */
void insert_robin_hood(RobinHoodTable *table, uint64_t kmer);

/**
 * @brief Counts k-mers in a given DNA sequence using the optimized Robin Hood
 * hash table.
 *
 * Iterates through the provided DNA sequence, extracts all k-mers of the
 * specified length 'k', and counts their occurrences, storing them in a Robin
 * Hood hash table.
 *
 * @param sequence The DNA sequence string to process.
 * @param k The length of the k-mers to count.
 * @return A pointer to a RobinHoodTable containing the k-mer counts, or NULL if
 * memory allocation fails.
 */
RobinHoodTable *cnt_kmer(const char *sequence, int k);

/**
 * @brief Retrieves the count of a specific k-mer from the Robin Hood hash
 * table.
 *
 * Looks up a given k-mer in the hash table and returns its associated count.
 * If the k-mer is not found, it returns 0.
 *
 * @param table A pointer to the RobinHoodTable to search.
 * @param kmer The k-mer (encoded as `uint64_t`) to look up.
 * @return The count of the k-mer if found, otherwise 0.
 */
uint32_t get_count_robin_hood(RobinHoodTable *table, uint64_t kmer);

/**
 * @brief Calculates the cosine similarity between two k-mer count tables.
 *
 * Computes the cosine similarity metric between two Robin Hood hash tables,
 * effectively comparing the k-mer frequency profiles of the two sequences
 * from which the tables were generated.
 *
 * @param table1 A pointer to the first RobinHoodTable.
 * @param table2 A pointer to the second RobinHoodTable.
 * @return The cosine similarity score (a double between 0.0 and 1.0).
 */
double cos_similarity(RobinHoodTable *table1,
                                   RobinHoodTable *table2);

/**
 * @brief Frees all memory allocated for a Robin Hood hash table.
 *
 * Releases the memory occupied by the hash table's data array and the table
 * structure itself. This function should be called when the table is no longer
 * needed to prevent memory leaks.
 *
 * @param table A pointer to the RobinHoodTable to be freed.
 */
void free_robin_hood_table(RobinHoodTable *table);

#endif // KMER_COUNTER_H
