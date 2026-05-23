#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>
#include <set>
#include <map>
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

int main() {

    /*Definiranje strukture ulaza (ulaz su .fastaq datoteke), granice za ulaz u klaster i liste klastera*/
    string header, NucleotideSequence, plus, qualityScore;
    int threshold = 12;
    vector<Cluster> myClusters;
    vector<Read> Reads; 

    /* Load all J*.fastq files */
    for (const auto& entry : fs::directory_iterator(".")) {

        string filename = entry.path().filename().string();

        /* Only process files starting with J */
        if (filename[0] == 'J') {

            cout << "Loading file: " << filename << endl;

            ifstream file(filename);

            while (
                getline(file, header) &&
                getline(file, NucleotideSequence) &&
                getline(file, plus) &&
                getline(file, qualityScore)
            ) {

                if (NucleotideSequence.length() == 296) {

                    Reads.push_back({
                        filename + " | " + header,
                        NucleotideSequence,
                        ""
                    });
                }
            }
        }
    }

    /* Generiranje MSA*/

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

    /* Klasteriranje */

    for (const Read& r : Reads) {

        if (myClusters.empty()) {
            Cluster newCluster;
            newCluster.sequences.push_back(r);
            myClusters.push_back(newCluster);
            continue;
        }

        int bestClusterIndex = -1;
        int bestDistance = INT_MAX;

        for (int i = 0; i < myClusters.size(); i++) {

            int clusterMin = INT_MAX;

            for (const Read& s : myClusters[i].sequences) {
                int dist = getHammingDistance(s.MSA_Sequence, r.MSA_Sequence);
                clusterMin = min(clusterMin, dist);
            }

            if (clusterMin < bestDistance) {
                bestDistance = clusterMin;
                bestClusterIndex = i;
            }
        }

        if (bestDistance <= threshold) {
            myClusters[bestClusterIndex].sequences.push_back(r);
        } else {
            Cluster newCluster;
            newCluster.sequences.push_back(r);
            myClusters.push_back(newCluster);
        }
    }

    for (Cluster& c : myClusters) {

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

    /*Ispis rezultata*/

    cout << "Clusters: " << myClusters.size() << endl;

    for (Cluster &c : myClusters){
        cout << "Cluster size: " << c.sequences.size() << endl;
        cout << "Representative: " << c.representativeSequence << endl << endl;
    }

    return 0;
}