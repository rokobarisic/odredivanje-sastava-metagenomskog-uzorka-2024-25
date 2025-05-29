#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta_parser.h"
#include "fastq_parser.h"
#include "kmer_counter.h"
// Compile flags: gcc -O3 -march=native -ffast-math -funroll-loops main.c -lm

int main(int argc, char **argv) {
  int kmer_length = 10;
  char *reference_filename = NULL;
  char *reads_filename = NULL;

  if (argc < 6) {
    fprintf(stderr,
            "Usage: %s -k <kmer_length> -ref <fasta_reference_file> -reads "
            "<fastq_reads_file>\n",
            argv[0]);
    fprintf(stderr,
            "  <fasta_reference_file>: Path to the FASTA reference genome.\n");
    fprintf(stderr, "  <fastq_reads_file>: Path to the FASTQ reads file.\n");
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-k")) {
      if (i + 1 < argc) {
        kmer_length = atoi(argv[i + 1]);
        if (kmer_length <= 0 || kmer_length > MAX_KMER_LEN) {
          fprintf(
              stderr,
              "Error: K-mer length must be a positive integer and max %d.\n",
              MAX_KMER_LEN);
          return 1;
        }
        i++;
      } else {
        fprintf(stderr, "Error: -k option requires an argument.\n");
        return 1;
      }
    } else if (!strcmp(argv[i], "-ref")) {
      if (i + 1 < argc) {
        reference_filename = argv[i + 1];
        i++;
      } else {
        fprintf(stderr, "Error: -ref option requires an argument.\n");
        return 1;
      }
    } else if (!strcmp(argv[i], "-reads")) {
      if (i + 1 < argc) {
        reads_filename = argv[i + 1];
        i++;
      } else {
        fprintf(stderr, "Error: -reads option requires an argument.\n");
        return 1;
      }
    } else {
      fprintf(stderr, "Error: Unknown option '%s'\n", argv[i]);
      fprintf(stderr,
              "Usage: %s -k <kmer_length> -ref <fasta_reference_file> -reads "
              "<fastq_reads_file>\n",
              argv[0]);
      return 1;
    }
  }

  if (reference_filename == NULL || reads_filename == NULL) {
    fprintf(stderr, "Error: Both reference file (-ref) and reads file (-reads) "
                    "are required.\n");
    fprintf(stderr,
            "Usage: %s -k <kmer_length> -ref <fasta_reference_file> -reads "
            "<fastq_reads_file>\n",
            argv[0]);
    return 1;
  }

  printf("K-mer length: %d\n", kmer_length);
  printf("Reference file: %s\n", reference_filename);
  printf("Reads file: %s\n", reads_filename);

  // --- PARSING FASTA REFERENCE FILE ---
  int num_ref_entries = 0;
  FastaEntry *ref_entries = parse_fasta(reference_filename, &num_ref_entries);
  if (ref_entries == NULL || num_ref_entries == 0) {
    fprintf(stderr,
            "Error: Failed to parse reference FASTA file or it's empty.\n");
    return 1;
  }
  printf("Successfully parsed %d FASTA reference entries.\n", num_ref_entries);

  // Use optimized Robin Hood hash table instead of KmerTable
  RobinHoodTable *reference_kmer_table = NULL;
  if (num_ref_entries > 0) {
    printf("Counting k-mers for reference genome: %s\n", ref_entries[0].id);
    reference_kmer_table =
        count_kmers_optimized(ref_entries[0].sequence, kmer_length);
    if (reference_kmer_table == NULL) {
      fprintf(stderr, "Error: Failed to count k-mers for reference.\n");
      free_fasta(ref_entries, num_ref_entries);
      return 1;
    }
    printf("Reference k-mer table created with %zu unique k-mers.\n",
           reference_kmer_table->size);
  } else {
    fprintf(stderr, "Error: No reference sequence found in FASTA file.\n");
    free_fasta(ref_entries, num_ref_entries);
    return 1;
  }

  // --- PARSING FASTQ READS FILE ---
  int num_read_entries = 0;
  char **read_entries = parse_fastq(reads_filename, &num_read_entries);
  if (read_entries == NULL || num_read_entries == 0) {
    fprintf(stderr, "Error: Failed to parse reads FASTQ file or it's empty.\n");
    free_fasta(ref_entries, num_ref_entries);
    free_robin_hood_table(reference_kmer_table);
    return 1;
  }
  printf("Successfully parsed %d FASTQ read entries.\n", num_read_entries);

  // --- OPTION 1: Process reads individually ---
  printf("\nCalculating cosine similarities for each read:\n");
  for (int i = 0; i < num_read_entries; i++) {
    RobinHoodTable *read_kmer_table =
        count_kmers_optimized(read_entries[i], kmer_length);
    double similarity =
        cosine_similarity_optimized(reference_kmer_table, read_kmer_table);
    printf("Read %d: Cosine Similarity = %.4f (k-mers: %zu)\n", i + 1,
           similarity, read_kmer_table->size);

    free_robin_hood_table(read_kmer_table);
  }

  // --- OPTION 2: Batch processing (uncomment to use instead of Option 1) ---
  // Implementation so far uses too much RAM and system kills the process.
  // Proper dynamic handling of RAM usage is needed to make such algorithm work.
  // printf("\nProcessing reads in batch mode...\n");
  // BatchResult* batch_result = process_reads_batch(reference_kmer_table,
  // read_entries, num_read_entries, kmer_length); if (batch_result) {
  //     printf("Batch processing complete. Results:\n");
  //     for (int i = 0; i < num_read_entries; i++) {
  //         printf("Read %d: Cosine Similarity = %.4f\n", i + 1,
  //         batch_result->similarities[i]);
  //     }
  //     free_batch_result(batch_result);
  // } else {
  //     fprintf(stderr, "Error: Batch processing failed.\n");
  // }

  // --- MEMORY CLEANUP ---
  free_fasta(ref_entries, num_ref_entries);
  free_robin_hood_table(reference_kmer_table);
  free_reads(read_entries, num_read_entries);

  printf("\nProgram finished successfully.\n");

  return 0;
}