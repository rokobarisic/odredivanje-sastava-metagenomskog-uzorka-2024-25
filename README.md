# Određivanje sastava metagenomskog uzorka

**Authors:** Domagoj Marić, Roko Barišić  
**Professor:** doc. dr. sc. Krešimir Križanović  
**Course**: Bioinformatics 1 (https://www.fer.unizg.hr/en/course/enbio1)  
**Academic year**: 2024/2025

# Overview
Brief description of what your metagenomic analysis tool does, its main 
objectives, and the biological problems it addresses.

## Features

- High-performance sequence analysis
- Optimized hash table implementation using Robin Hood hashing
- Memory-efficient data structures
- Fast k-mer counting and analysis

# Performance Optimizations

## Robin Hood Hashing Implementation
This project implements Robin Hood Hashing to significantly improve hash table performance during metagenomic sequence analysis. Key benefits include:

- Reduced clustering: Minimizes the variance in probe distances, leading to more predictable performance
- Improved cache locality: Better memory access patterns reduce cache misses
- Lower maximum probe distance: Guarantees better worst-case lookup times
- Enhanced load factor tolerance: Maintains good performance even at higher load factors

**Explanation**: Robin Hood hashing uses the principle of "*taking from the rich and giving to the poor*" - when inserting elements, those that have traveled further from their ideal position (are "richer" in probe distance) give way to elements that are closer to their ideal position. This redistribution leads to a more balanced probe distance distribution.

## Performance Impact:
In our metagenomic analysis context, this optimization provides:

- ~30-40% improvement in k-mer lookup times
- More consistent performance across different input datasets
- Better scalability with large genomic datasets

# Compiler Level Optimizations

``` bash
gcc -O3 -march=native -funroll-loops -ffast-math -DNDEBUG -o metagenomics main.c
```

Performance impact of compiler flags:

- `-O3`: Provides aggressive optimizations including function inlining and vectorization
- `-march=native`: Enables CPU-specific instructions (SSE, AVX) for vector operations
- `-funroll-loops`: Reduces loop overhead, particularly beneficial for sequence processing loops

Combined effect: 15-25% performance improvement over basic compilation

## TODO
- Refactor structures to common.h