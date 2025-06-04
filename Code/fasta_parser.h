#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @def MAX_SEQ
 * @brief Maximum expected line length in a FASTA file.
 *
 * Defines the maximum buffer size for reading a single line
 * from a FASTA file. Adjust this value as needed.
 */
#define MAX_SEQ 1024

/**
 * @struct FastaEntry
 * @brief Stores a single FASTA sequence entry.
 *
 * Contains the ID (header) and the sequence for an individual
 * entry in a FASTA file. All strings are dynamically allocated
 * and must be freed by the user.
 */
typedef struct {
  char *id;       /**< @brief Sequence ID(e.g., ">accession_id description"). */
  char *sequence; /**< @brief The nucleotide or protein sequence. */
} FastaEntry;

/**
 * @brief Parses FASTA file and loads entries into memory.
 *
 * Reads a FASTA formatted file, allocates memory for each sequence and its ID
 * and stores them in an array of FastaEntry structures. User is responsible
 * for freeing the allocated memory by calling the `free_refs` function.
 *
 * @param filename The path to the FASTA file.
 * @param count A pointer to an integer where the number of parsed entries will
 * be stored.
 * @return A pointer to a dynamically allocated array of FastaEntry structures,
 * or NULL in case of an error (e.g., file cannot be opened).
 */
FastaEntry *parse_fasta(const char *filename, int *count);

/**
 * @brief Frees memory allocated for FASTA entries.
 *
 * Iterates through an array of FastaEntry structures and frees the memory
 * allocated for the IDs and sequences within each structure,
 * and finally frees the array of structures itself.
 *
 * @param entries A pointer to the array of FastaEntry structures to be freed.
 * @param count The number of entries in the 'entries' array.
 */
void free_refs(FastaEntry *entries, int count);

#endif // FASTA_PARSER_H
