#include <dirent.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "fasta_parser.h"
#include "fastq_parser.h"
#include "kmer_counter.h"

#define MAX_REF_FILES 64

int has_extension(const char *filename, const char *ext) {
  const char *dot = strrchr(filename, '.');
  return dot && strcmp(dot, ext) == 0;
}

typedef struct {
  char name[1024];
  RobinHoodTable *table;
} ReferenceData;

void process_files(const char *ref_dir, const char *reads_dir, int kmer_length,
                   double threshold) {
  DIR *ref_dp = opendir(ref_dir);
  if (!ref_dp) {
    perror("Error opening reference directory");
    exit(1);
  }

  // Load all reference files
  ReferenceData references[MAX_REF_FILES];
  int num_refs = 0;

  struct dirent *ref_entry;
  while ((ref_entry = readdir(ref_dp)) != NULL) {
    if (ref_entry->d_type != DT_REG ||
        !(has_extension(ref_entry->d_name, ".fasta") ||
          has_extension(ref_entry->d_name, ".fa")))
      continue;

    char ref_path[1024];
    snprintf(ref_path, sizeof(ref_path), "%s/%s", ref_dir, ref_entry->d_name);

    int num_entries = 0;
    FastaEntry *entries = parse_fasta(ref_path, &num_entries);
    if (!entries || num_entries == 0) {
      fprintf(stderr, "Warning: Skipping empty or invalid FASTA: %s\n",
              ref_path);
      continue;
    }

    RobinHoodTable *kmer_table = cnt_kmer(entries[0].sequence, kmer_length);
    if (!kmer_table) {
      fprintf(stderr, "Error: k-mer counting failed for %s\n", ref_path);
      free_refs(entries, num_entries);
      continue;
    }

    strncpy(references[num_refs].name, ref_entry->d_name,
            sizeof(references[num_refs].name));
    references[num_refs].table = kmer_table;
    num_refs++;

    free_refs(entries, num_entries);
    if (num_refs >= MAX_REF_FILES) {
      fprintf(stderr, "Warning: Max reference limit reached.\n");
      break;
    }
  }
  closedir(ref_dp);

  // Process reads
  DIR *reads_dp = opendir(reads_dir);
  if (!reads_dp) {
    perror("Error opening reads directory");
    exit(1);
  }

  struct dirent *reads_entry;
  while ((reads_entry = readdir(reads_dp)) != NULL) {
    if (reads_entry->d_type != DT_REG ||
        !(has_extension(reads_entry->d_name, ".fastq") ||
          has_extension(reads_entry->d_name, ".fq")))
      continue;

    char reads_path[1024];
    snprintf(reads_path, sizeof(reads_path), "%s/%s", reads_dir,
             reads_entry->d_name);

    int num_reads = 0;
    char **read_entries = parse_fastq(reads_path, &num_reads);
    if (!read_entries || num_reads == 0) {
      fprintf(stderr, "Warning: Skipping empty/invalid FASTQ: %s\n",
              reads_path);
      continue;
    }

    int *best_ref_idx = calloc(num_reads, sizeof(int));
    double *best_sim_val = calloc(num_reads, sizeof(double));
    int passing_reads = 0;

    for (int i = 0; i < num_reads; i++) {
      best_ref_idx[i] = -1;
      best_sim_val[i] = 0.0;

      RobinHoodTable *read_table = cnt_kmer(read_entries[i], kmer_length);
      if (!read_table)
        continue;

      for (int r = 0; r < num_refs; r++) {
        double sim = cos_similarity(references[r].table, read_table);
        if (sim > best_sim_val[i]) {
          best_sim_val[i] = sim;
          best_ref_idx[i] = r;
        }
      }

      if (best_sim_val[i] >= threshold)
        passing_reads++;

      free_robin_hood_table(read_table);
    }

    float portion = (1.0 * passing_reads) / num_reads;
    // Summary Report
    printf("\n=== Summary for Reads: %s ===\n", reads_entry->d_name);
    printf("Total reads: %d\n", num_reads);
    printf("Reads above threshold: %d (%.4f)\n", passing_reads, portion);

    free(best_ref_idx);
    free(best_sim_val);
    free_reads(read_entries, num_reads);
  }

  closedir(reads_dp);

  for (int r = 0; r < num_refs; r++)
    free_robin_hood_table(references[r].table);
}

int main(int argc, char **argv) {
  int kmer_length = 5;
  char *ref_dir = NULL;
  char *reads_dir = NULL;
  double threshold = 0.0;

  if (argc < 6) {
    fprintf(stderr,
            "Usage: %s -k <kmer_length> -refdir <fasta_dir> -readsdir "
            "<fastq_dir> [-t <threshold>]\n",
            argv[0]);
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-k") && i + 1 < argc)
      kmer_length = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-refdir") && i + 1 < argc)
      ref_dir = argv[++i];
    else if (!strcmp(argv[i], "-readsdir") && i + 1 < argc)
      reads_dir = argv[++i];
    else if (!strcmp(argv[i], "-t") && i + 1 < argc)
      threshold = atof(argv[++i]);
    else {
      fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
      return 1;
    }
  }

  if (!ref_dir || !reads_dir || kmer_length <= 0 ||
      kmer_length > MAX_KMER_LEN) {
    fprintf(stderr, "Error: Missing or invalid arguments.\n");
    return 1;
  }

  printf(
      "K-mer length: %d\nReference dir: %s\nReads dir: %s\nThreshold: %.3f\n",
      kmer_length, ref_dir, reads_dir, threshold);

  process_files(ref_dir, reads_dir, kmer_length, threshold);

  return 0;
}
