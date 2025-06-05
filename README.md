# Determination of the composition of a metagenomic sample

**Authors:** Domagoj Marić, Roko Barišić  
**Professor:** doc. dr. sc. Krešimir Križanović  
**Course**: Bioinformatics 1 (https://www.fer.unizg.hr/en/course/enbio1)  
**Academic year**: 2024/2025

## Features

- 🚀 High-performance sequence analysis optimized for large-scale FASTA and FASTQ files
- 🧠 Memory-efficient data structures for k-mer counting
- ⚙️ Robin Hood Hashing for fast and balanced hash table operations
- 🔍 Cosine similarity between k-mer vectors for comparative genomics
- 🧬 Supports arbitrary k-mer lengths up to 31
- 📚 No external libraries required


# Performance Optimizations

We decided to use **Robin Hood Hashing** to optimize hash table behavior in 
high-throughput genomic contexts.

## Benefits

- 📉 Reduced clustering → faster lookups
- 🧠 Improved cache locality → lower memory latency
- ⏱️ Lower worst-case probe time
- 📈 High load factor tolerance → handles larger tables before resizing

**Explanation**: Robin Hood hashing uses the principle of "*taking from the rich and giving to the poor*" - when inserting elements, those that have traveled further from their ideal position (are "richer" in probe distance) give way to elements that are closer to their ideal position. This redistribution leads to a more balanced probe distance distribution.

## Performance Impact:
In our metagenomic analysis context, this optimization provides:

- ~30-40% improvement in k-mer lookup times
- Better handling of high-throughput reads (e.g. >500,000 FASTQ reads)
- Better scalability with large genomic datasets

## Compiler Level Optimizations

``` bash
gcc -O3 -march=native -funroll-loops -ffast-math *.c -o metagenomics
```

Performance impact of compiler flags:

- `-O3`: Provides aggressive optimizations including function inlining and vectorization
- `-march=native`: Enables CPU supported instructions (SSE, AVX) for vector operations
- `-funroll-loops`: Reduces loop overhead, particularly beneficial for sequence processing loops
- `-ffast-math`: Enables faster (non-strict ones) floating-point math

Combined effect: 15-25% performance improvement over basic compilation

# Example usage

Precompiled binaries are available on our GitHub Releases page.
The program accepts directories containing FASTA reference genomes and FASTQ read files
with *optional* threshold parameter defining similarity condition, e.g. threshold  > 0.5:

```bash
./metagenomics -k <kmer_length> -refdir <fasta_reference_dir> -readsdir <fastq_reads_dir> [-t <threshold>]
```

In a benchmark using 4 reference genomes and 4 read files (~6,2GB total size), the
program completed in under 2 minutes on a machine with Intel Core i5 6500. Benchmark output:

```
K-mer length: 5
Reference dir: ../Data/References/
Reads dir: ../Data/Readings/
Threshold: 0.100

=== Summary for Reads: NCTC12874.fastq ===
Total reads: 344672
Reads above threshold: 344154 (0.9985)

=== Summary for Reads: NCTC10116.fastq ===
Total reads: 561896
Reads above threshold: 561546 (0.9994)

=== Summary for Reads: NCTC11829.fastq ===
Total reads: 118225
Reads above threshold: 118165 (0.9995)

=== Summary for Reads: NCTC12157.fastq ===
Total reads: 209309
Reads above threshold: 208860 (0.9979)
```

## Notes on K-mer Length and Similarity

- Not all possible k-mers are observed in real genomes. For example, with `k = 7`, there are 16,384 possible
combinations (4^7), but actual k-mer counts are often slightly lower due to natural sequence composition.
- Shorter k-mers (e.g., `k = 5`) tend to result in higher cosine similarity scores due to higher chance overlap.
- Longer k-mers are more specific but rare. They require near-perfect alignment, so similarities drop — this
is expected and biologically meaningful.

## TODO
- Utilize multiprocessing for comparisons
- Export results in standard formats (CSV/JSON)

## License
This project is open source under the MIT License.
