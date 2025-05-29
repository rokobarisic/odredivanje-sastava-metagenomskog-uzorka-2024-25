#ifndef FASTQ_PARSER_H
#define FASTQ_PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char **parse_fastq(const char *filename, int *count);
void free_reads(char **reads, int count);

inline char **parse_fastq(const char *filename, int *count) {
  FILE *file = fopen(filename, "r");
  if (!file)
    return NULL;

  int capacity = 100;
  int index = 0;
  char **reads = (char **)malloc(capacity * sizeof(char *));

  char line[1024];
  int line_num = 0;
  while (fgets(line, sizeof(line), file)) {
    line_num++;
    if (line_num % 4 == 2) {
      line[strcspn(line, "\n")] = 0;
      reads[index] = strdup(line);
      index++;
      if (index >= capacity) {
        capacity *= 2;
        reads = (char **)realloc(reads, capacity * sizeof(char *));
      }
    }
  }

  fclose(file);
  *count = index;
  return reads;
}

inline void free_reads(char **reads, int count) {
  for (int i = 0; i < count; i++) {
    free(reads[i]);
  }
  free(reads);
}

#endif // FASTQ_PARSER_H
