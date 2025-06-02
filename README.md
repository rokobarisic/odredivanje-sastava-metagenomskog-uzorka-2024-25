# Determination of the composition of a metagenomic sample

**Authors:** Domagoj Marić, Roko Barišić  
**Professor:** doc. dr. sc. Krešimir Križanović  
**Course**: Bioinformatics 1 (https://www.fer.unizg.hr/en/course/enbio1)  
**Academic year**: 2024/2025

# How to use?

1. git clone https://github.com/rokobarisic/odredivanje-sastava-metagenomskog-uzorka-2024-25.git
2. go to the folder where you cloned the repository and position yourself in its subfolder "external"
3. git clone https://github.com/rvaser/bioparser
4. put reference genome files in /Data/References
5. put readings in /Data/Readings
6. go to folder Code
7. g++ reading_creator.cpp -o reading_creator
8. ./reading_creator
9. g++ main.cpp -std=c++17 -I../external/bioparser/include -o main -lz
10. ./main
