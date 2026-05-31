# bioinformatics-allele-clustering
Allele clustering project for Bioinformatics 1

# How to use

 1. Clone repository\n
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name

 2. Build the project with: 
g++ Bioinformatika.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o Bioinformatika

 3. Run the program
./Bioinformatika

# Input data
The program expects either a single .fastaq file or a folder path containing FASTQ files.

Example structure:
fastq/
 ├── J1_S_CE_IonXpress_018.fastq
 ├── J2_S_CE_IonXpress_019.fastq
 ├── ...

 # Requirements
- C++17 or higher
- g++ compiler
- filesystem support (std::filesystem)
- SPOA library (https://github.com/rvaser/spoa)
