#include "fastq_parser.h"

char **parse_fastq(const char *filename, int *count) {
  FILE *file = fopen(filename, "r");
  if (!file)
    return NULL;

  int capacity = 100;
  int index = 0;
  char **reads = (char **)malloc(capacity * sizeof(char *));

  char line[1024];
  int line_num = 0;

  while (fgets(line, sizeof(line), file)) {
    if (line[0] != '@')
      continue; // Start of a read

    if (!fgets(line, sizeof(line), file))
      break; // Sequence line
    line[strcspn(line, "\n")] = 0;

    reads[index] = strdup(line);
    index++;
    if (index >= capacity) {
      capacity *= 2;
      reads = (char **)realloc(reads, capacity * sizeof(char *));
    }

    // Skip '+' line and quality line
    fgets(line, sizeof(line), file); // '+' line
    fgets(line, sizeof(line), file); // quality line

    for (char *p = reads[index - 1]; *p; ++p) {
      *p = toupper(*p);
      if (*p != 'A' && *p != 'C' && *p != 'G' && *p != 'T') {
        *p = 'N'; // or remove/skip this read later
      }
    }
  }

  fclose(file);
  *count = index;
  return reads;
}

void free_reads(char **reads, int count) {
  for (int i = 0; i < count; i++) {
    free(reads[i]);
  }
  free(reads);
}