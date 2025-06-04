#ifndef FASTQ_PARSER_H
#define FASTQ_PARSER_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Parses FASTQ file and extracts sequence reads.
 *
 * Reads a FASTQ formatted file, extracts each sequence read (the second line of
 * each FASTQ entry), and stores them as dynamically allocated strings in a
 * dynamically allocated array of char pointers. User is responsible for freeing
 * the allocated memory by calling the `free_reads` function. Note: This parser
 * currently only extracts the sequence string and does not handle IDs or
 * quality scores.
 *
 * @param filename The path to the FASTQ file.
 * @param count A pointer to an integer where the total number of parsed reads
 * will be stored.
 * @return A pointer to a dynamically allocated array of char pointers, where
 * each pointer points to a dynamically allocated string representing a sequence
 * read. Returns NULL if the file cannot be opened or if memory allocation
 * fails.
 */
char **parse_fastq(const char *filename, int *count);

/**
 * @brief Frees allocated memory for FASTQ reads.
 *
 * This function iterates through the array of char pointers (reads) and frees
 * the memory allocated for each individual sequence string, and then frees
 * the array of pointers itself. This should be called after `parse_fastq`
 * when the reads are no longer needed to prevent memory leaks.
 *
 * @param reads A pointer to the array of char pointers (sequence reads) to be
 * freed.
 * @param count The number of reads in the 'reads' array. This should be the
 * value returned by the 'count' parameter of `parse_fastq`.
 */
void free_reads(char **reads, int count);

#endif // FASTQ_PARSER_H
