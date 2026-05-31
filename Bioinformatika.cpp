#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>
#include <set>
#include <map>
#include <algorithm>

#include "spoa/include/spoa/spoa.hpp"

#include <filesystem>
namespace fs = std::filesystem;

using namespace std;

/* Očitanja s nazivom, originalnom sekvencom i poravnatom sekvencom */
struct Read{
    string name;
    string originalSequence;
    string MSA_Sequence;
};

/*Struktura podataka Cluster za spremanje sličnih sekvenci, sastoji se od vektora sekvenci i jedne konsenzusne sekvence*/
struct Cluster{
    string representativeSequence;
    vector<Read> sequences;
};

/* Funkcija za računanje Hammingove udaljenosti između dvije sekvence*/
int getHammingDistance(const string& sequence1, const string& sequence2){
    int distance = 0;
    for (size_t i = 0; i < sequence1.size(); i++) {
        if (sequence1[i] != sequence2[i]) distance++;
    }
    return distance;
}

/*Generiranje MSA*/
void GenerateMSA(vector<Read> &Reads){
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 0, -1, -1);  

    spoa::Graph graph{};

    for (const auto& it : Reads) {
        auto alignment = alignment_engine->Align(it.originalSequence, graph);
        graph.AddAlignment(alignment, it.originalSequence);
    }

    auto msa = graph.GenerateMultipleSequenceAlignment();

    for (size_t i = 0; i < Reads.size(); i++) {
        Reads[i].MSA_Sequence = msa[i];
    }

    return;
}

/*Klasteriranje - prolazimo po klasterima i svim njihovim trenutnim članovima (sekvencama) kako bi našli onu s najmanjom udaljnosti i zatim 
sekvencu dodajemo u klaster u kojem se nalazi sekvenca s najmanjom razlikom u udaljenosti s time da ne smije prelaziti threshold*/
vector<Cluster> Clustering(vector<Read> &Reads, int threshold){

    vector<Cluster> Clusters;

    for (const Read& r : Reads) {

        if (Clusters.empty()) {
            Cluster newCluster;
            newCluster.sequences.push_back(r);
            Clusters.push_back(newCluster);
            continue;
        }

        int bestClusterIndex = -1;
        int bestDistance = INT_MAX;

        for (int i = 0; i < Clusters.size(); i++) {

            int clusterMin = INT_MAX;

            for (const Read& s : Clusters[i].sequences) {
                int dist = getHammingDistance(s.MSA_Sequence, r.MSA_Sequence);
                clusterMin = min(clusterMin, dist);
            }

            if (clusterMin < bestDistance) {
                bestDistance = clusterMin;
                bestClusterIndex = i;
            }
        }

        if (bestDistance <= threshold) {
            Clusters[bestClusterIndex].sequences.push_back(r);
        } else {
            Cluster newCluster;
            newCluster.sequences.push_back(r);
            Clusters.push_back(newCluster);
        }
    }
    
    return Clusters;
}

/*Generiranje konsenzusne sekvence*/
void GenerateConsensus(vector<Cluster> &clusters){
    for (Cluster& c : clusters) {

        auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 0, -1, -1
        );

        spoa::Graph graph;

        for (const Read& r : c.sequences) {
            auto alignment = alignment_engine->Align(r.originalSequence, graph);
            graph.AddAlignment(alignment, r.originalSequence);
        }

        auto consensus = graph.GenerateConsensus();
        c.representativeSequence = consensus;
    }
}

/*Učitavanje podataka*/
vector<Read> LoadReads(const fs::path& filepath){
    vector<Read> reads;

    string header, sequence, plus, quality;

    ifstream file(filepath);

    /*Učitavamo po četiri retka odjednom naziv očitanja, sekvencu, znak plus i 'quality score' iz .fastaq datoteke*/
    while (
        getline(file, header) &&
        getline(file, sequence) &&
        getline(file, plus) &&
        getline(file, quality)
    ) {
        if (sequence.length() == 296) {
            reads.push_back({
                filepath.filename().string() + " | " + header,
                sequence,
                ""
            });
        }
    }

    return reads;
}

/* Glavni program */

int main() {

    string inputPath;

    cout << "Enter FASTQ file or folder: ";
    getline(cin, inputPath);

    bool singleFile = false;

    if (fs::is_regular_file(inputPath)) {
        singleFile = true;
    } else if (!fs::is_directory(inputPath)) {
        cerr << "Invalid path!" << endl;
        return 1;
    }

    /*Definiranje granice za ulaz u klaster*/
    int threshold = 12;

    /*Način rada za samo jednu datoteku*/
    if (singleFile){
        
        ifstream file(inputPath);

        vector<Read> Reads = LoadReads(inputPath);

        GenerateMSA(Reads);
        vector<Cluster> myClusters = Clustering(Reads, threshold);
        
        /*Izbacujemo male klastere*/
        myClusters.erase(
            remove_if(
                    myClusters.begin(),
                    myClusters.end(),
                        [](const Cluster& c) {
                        return c.sequences.size() <= 2;
                    }
                ),
                    myClusters.end()
        );

        GenerateConsensus(myClusters);    

        /*Ispis rezultata u datoteku u FASTA formatu i membership.txt*/
        ofstream fasta(fs::path(inputPath).stem().string() + "_variants.fasta");

        for (size_t i = 0; i < myClusters.size(); i++) {

            fasta << ">Variant_" << i + 1
                << "_size_" << myClusters[i].sequences.size()
                << "\n";

            fasta << myClusters[i].representativeSequence
                << "\n";
        }

        ofstream clusters(fs::path(inputPath).stem().string() + "_cluster_membership.txt");

        for (size_t i = 0; i < myClusters.size(); i++) {

            clusters << "Cluster "
                    << i + 1
                    << "\n";

            clusters << "Representative:\n";
            clusters << myClusters[i].representativeSequence
                    << "\n\n";

            clusters << "Reads:\n";

            for (const Read& r : myClusters[i].sequences) {
                clusters << r.name << "\n";
            }

            clusters << "\n-----------------------------------\n\n";
        }

    /*Način rada da folder s datotekama*/
    } else {

        vector<Cluster> allClusters;
        vector<Read> allRepresentatives;

        for (const auto& entry : fs::directory_iterator(inputPath)){

            cout << "Loading file: " << entry.path().filename().string() << endl;

            vector<Read> Reads = LoadReads(entry.path());

            GenerateMSA(Reads);
            vector<Cluster> myClusters = Clustering(Reads, threshold);

            myClusters.erase(
                remove_if(
                    myClusters.begin(),
                    myClusters.end(),
                        [](const Cluster& c) {
                        return c.sequences.size() <= 2;
                    }
                ),
                    myClusters.end()
            );

            GenerateConsensus(myClusters);

            for (Cluster& c : myClusters) {
                Read newRepresentative;
                newRepresentative.name = "Variant_representative_of_" + entry.path().filename().string();
                newRepresentative.originalSequence = c.representativeSequence;
                allRepresentatives.push_back(newRepresentative);
            }
        }

        /*Klasteriranje dobivenih varijanti iz svakog uzorka*/
        GenerateMSA(allRepresentatives);

        allClusters = Clustering(allRepresentatives, threshold);
        GenerateConsensus(allClusters);

        ofstream fasta("Variants.fasta");

        for (size_t i = 0; i < allClusters.size(); i++) {

            fasta << ">Variant_" << i + 1
                << "_size_" << allClusters[i].sequences.size()
                << "\n";

            fasta << allClusters[i].representativeSequence
                << "\n";
        }

        ofstream clusters("Cluster_membership.txt");

        for (size_t i = 0; i < allClusters.size(); i++) {

            clusters << "Cluster "
                    << i + 1
                    << "\n";

            clusters << "Representative:\n";
            clusters << allClusters[i].representativeSequence
                    << "\n\n";

            clusters << "Reads:\n";

            for (const Read& r : allClusters[i].sequences) {
                clusters << r.name << "\n";
            }

            clusters << "\n--------------------------------\n\n";
        }

    }

    return 0;
}