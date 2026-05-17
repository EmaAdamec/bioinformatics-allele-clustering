#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "spoa/include/spoa/spoa.hpp"

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

    /*Učitavanje datoteke*/
    ifstream file("J29_B_CE_IonXpress_005.fastq");

    /*Čitanje datoteke, po 4 reda odjednom: header, sekvenca, plus i ocjena kvalitete (gleda se samo sekvenca, ostala tri reda zanemarujemo)*/
    while (getline(file, header) && getline(file, NucleotideSequence) && getline(file, plus) && getline(file, qualityScore)){

        /*Gledamo samo sekvence najčešće duljine za sad*/
        if (NucleotideSequence.length() == 296) {
            Reads.push_back({
                header,
                NucleotideSequence,
                ""
            });
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

        bool assigned = false;

        for (Cluster& c : myClusters) {

            bool fitsCluster = true;

            for (const Read& s : c.sequences) {
                if (getHammingDistance(s.MSA_Sequence, r.MSA_Sequence) >= threshold) {
                    fitsCluster = false;
                    break;
                }
            }

            if (fitsCluster) {
                c.sequences.push_back(r);
                assigned = true;
                break;
            }
        }

        if (!assigned) {
            Cluster newCluster;
            newCluster.representativeSequence = "";
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