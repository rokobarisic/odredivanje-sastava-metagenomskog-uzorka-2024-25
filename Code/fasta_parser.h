#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ 1024

typedef struct {
  char *id;
  char *sequence;
} FastaEntry;

FastaEntry *parse_fasta(const char *filename, int *count);
void free_fasta(FastaEntry *entries, int count);

inline FastaEntry *parse_fasta(const char *filename, int *count) {
  FILE *file = fopen(filename, "r");
  if (!file)
    return NULL;

  int capacity = 10;
  int index = 0;
  FastaEntry *entries = (FastaEntry *)malloc(capacity * sizeof(FastaEntry));

  char line[MAX_SEQ];
  char *current_seq = NULL;
  char *current_id = NULL;
  size_t seq_len = 0;

  while (fgets(line, sizeof(line), file)) {
    if (line[0] == '>') {
      if (current_id) {
        entries[index].id = strdup(current_id);
        entries[index].sequence = strdup(current_seq);
        index++;
        if (index >= capacity) {
          capacity *= 2;
          entries =
              (FastaEntry *)realloc(entries, capacity * sizeof(FastaEntry));
        }
        free(current_seq);
        seq_len = 0;
      }
      line[strcspn(line, "\n")] = 0;
      current_id = line + 1;
      current_seq = (char *)calloc(1, sizeof(char));
    } else {
      line[strcspn(line, "\n")] = 0;
      size_t len = strlen(line);
      current_seq = (char *)realloc(current_seq, seq_len + len + 1);
      strcpy(current_seq + seq_len, line);
      seq_len += len;
    }
  }
  if (current_id && current_seq) {
    entries[index].id = strdup(current_id);
    entries[index].sequence = strdup(current_seq);
    index++;
    free(current_seq);
  }

  fclose(file);
  *count = index;
  return entries;
}

inline void free_fasta(FastaEntry *entries, int count) {
  for (int i = 0; i < count; i++) {
    free(entries[i].id);
    free(entries[i].sequence);
  }
  free(entries);
}

#endif // FASTA_PARSER_H