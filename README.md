# bioinformatics-allele-clustering
Allele clustering project for Bioinformatics 1

# How to use

# 1. Clone repository
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name

# 2. Build the project
g++ Bioinformatika.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o Bioinformatika

# 3. Run the program
./Bioinformatika

# Input data
The program expects either a single .fastaq file or a folder path containing FASTQ files.

 # Requirements
- C++17 or higher
- g++ compiler
- filesystem support (std::filesystem)
- SPOA library (https://github.com/rvaser/spoa)
