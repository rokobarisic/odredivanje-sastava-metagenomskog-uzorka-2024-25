#include "fasta_parser.h"

FastaEntry *parse_fasta(const char *filename, int *count) {
  FILE *file = fopen(filename, "r");
  if (!file)
    return NULL;

  int capacity = 10;
  int index = 0;
  FastaEntry *entries = (FastaEntry *)malloc(capacity * sizeof(FastaEntry));

  char line[MAX_SEQ];
  char *current_seq = NULL;
  char *current_id = NULL;
  size_t seq_capacity = 0;
  size_t seq_len = 0;

  while (fgets(line, sizeof(line), file)) {
    if (line[0] == '>') {
      // Store previous entry if it exists
      if (current_id) {
        entries[index].id = current_id;
        entries[index].sequence = current_seq;
        index++;
        if (index >= capacity) {
          capacity *= 2;
          entries = (FastaEntry *)realloc(entries, capacity * sizeof(FastaEntry));
        }
        current_id = NULL;
        current_seq = NULL;
        seq_len = 0;
        seq_capacity = 0;
      }

      line[strcspn(line, "\n")] = 0;
      current_id = strdup(line + 1);
    } else {
      line[strcspn(line, "\n")] = 0;
      size_t len = strlen(line);
      if (seq_len + len + 1 > seq_capacity) {
        seq_capacity = (seq_len + len + 1) * 2;
        current_seq = (char *)realloc(current_seq, seq_capacity);
      }
      memcpy(current_seq + seq_len, line, len);
      seq_len += len;
      current_seq[seq_len] = '\0';
    }
  }

  // Final entry
  if (current_id && current_seq) {
    entries[index].id = current_id;
    entries[index].sequence = current_seq;
    index++;
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