#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "fasta_parser.h"
#include "fastq_parser.h"
#include "kmer_counter.h"

int has_extension(const char *filename, const char *ext) {
  const char *dot = strrchr(filename, '.');
  return dot && strcmp(dot, ext) == 0;
}

void process_files(const char *ref_dir, const char *reads_dir,
                   int kmer_length) {
  DIR *ref_dp = opendir(ref_dir);
  DIR *reads_dp = opendir(reads_dir);

  if (!ref_dp) {
    perror("Error opening reference directory");
    exit(1);
  }
  if (!reads_dp) {
    perror("Error opening reads directory");
    exit(1);
  }

  struct dirent *ref_entry;
  while ((ref_entry = readdir(ref_dp)) != NULL) {
    if (ref_entry->d_type != DT_REG ||
        !(has_extension(ref_entry->d_name, ".fasta") ||
          has_extension(ref_entry->d_name, ".fa")))
      continue;

    char ref_path[1024];
    snprintf(ref_path, sizeof(ref_path), "%s/%s", ref_dir, ref_entry->d_name);

    int num_ref_entries = 0;

    FastaEntry *ref_entries = parse_fasta(ref_path, &num_ref_entries);
    if (!ref_entries || num_ref_entries == 0) {
      fprintf(stderr, "Warning: Skipping empty/invalid FASTA file: %s\n",
              ref_path);
      continue;
    }

    RobinHoodTable *reference_kmer_table =
        cnt_kmer(ref_entries[0].sequence, kmer_length);
    if (!reference_kmer_table) {
      fprintf(stderr, "Error: Failed to count k-mers for %s\n", ref_path);
      free_refs(ref_entries, num_ref_entries);
      continue;
    } else
      printf("Reference k-mer table created with %zu unique k-mers.\n",
             reference_kmer_table->size);

    struct dirent *reads_entry;
    rewinddir(reads_dp); // Restart reads directory iteration for each reference

    while ((reads_entry = readdir(reads_dp)) != NULL) {
      if (reads_entry->d_type != DT_REG ||
          !(has_extension(reads_entry->d_name, ".fastq") ||
            has_extension(reads_entry->d_name, ".fq")))
        continue;

      char reads_path[1024];
      snprintf(reads_path, sizeof(reads_path), "%s/%s", reads_dir,
               reads_entry->d_name);

      int num_read_entries = 0;
      char **read_entries = parse_fastq(reads_path, &num_read_entries);
      if (!read_entries || num_read_entries == 0) {
        fprintf(stderr, "Warning: Skipping empty/invalid FASTQ file: %s\n",
                reads_path);
        continue;
      }

      double max_similarity = 0.0;
      int max_idx = 0;

      for (int i = 0; i < num_read_entries; i++) {
        RobinHoodTable *read_kmer_table =
            cnt_kmer(read_entries[i], kmer_length);
        double similarity =
            cos_similarity(reference_kmer_table, read_kmer_table);
        if (similarity > max_similarity) {
          max_similarity = similarity;
          max_idx = i;
        }
        free_robin_hood_table(read_kmer_table);
      }

      printf("Reference: %s | Reads: %s | Max similarity: %.6lf (read #%d)\n",
             ref_entry->d_name, reads_entry->d_name, max_similarity, max_idx);

      free_reads(read_entries, num_read_entries);
    }

    free_refs(ref_entries, num_ref_entries);
    free_robin_hood_table(reference_kmer_table);
  }

  closedir(ref_dp);
  closedir(reads_dp);
}

int main(int argc, char **argv) {
  int kmer_length = 5;
  char *ref_dir = NULL;
  char *reads_dir = NULL;

  if (argc < 6) {
    fprintf(stderr,
            "Usage: %s -k <kmer_length> -refdir <fasta_directory> -readsdir "
            "<fastq_directory>\n",
            argv[0]);
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-k")) {
      if (i + 1 < argc) {
        kmer_length = atoi(argv[++i]);
        if (kmer_length <= 0 || kmer_length > MAX_KMER_LEN) {
          fprintf(stderr, "Error: Invalid k-mer length.\n");
          return 1;
        }
      } else {
        fprintf(stderr, "Error: -k requires an argument.\n");
        return 1;
      }
    } else if (!strcmp(argv[i], "-refdir")) {
      if (i + 1 < argc)
        ref_dir = argv[++i];
      else {
        fprintf(stderr, "Error: -refdir requires an argument.\n");
        return 1;
      }
    } else if (!strcmp(argv[i], "-readsdir")) {
      if (i + 1 < argc)
        reads_dir = argv[++i];
      else {
        fprintf(stderr, "Error: -readsdir requires an argument.\n");
        return 1;
      }
    } else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
    }
  }

  if (!ref_dir || !reads_dir) {
    fprintf(stderr, "Both -refdir and -readsdir must be specified.\n");
    return 1;
  }

  printf("K-mer length: %d\nReference dir: %s\nReads dir: %s\n\n", kmer_length,
         ref_dir, reads_dir);

  process_files(ref_dir, reads_dir, kmer_length);

  printf("\nProgram finished successfully.\n");

  return 0;
}