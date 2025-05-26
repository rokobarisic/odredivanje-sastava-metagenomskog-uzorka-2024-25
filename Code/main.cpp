#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include "bioparser/fasta_parser.hpp"
#include <math.h>

using namespace std;
namespace fs = std::filesystem;

struct Sequence
{
public:
    string id;
    string data;

    Sequence(const char *id_ptr, uint32_t id_len,
             const char *data_ptr, uint32_t data_len)
    {
        id = string(id_ptr, id_len);
        data = string(data_ptr, data_len);
    }

    void Print() const
    {
        cout << "ID: " << id << "\n";
        cout << "Sekvenca: " << data << "\n\n";
    }
};

uint16_t encode(const std::string &kmer)
{
    if (kmer.size() != 5)
        throw invalid_argument("Sekvenca mora biti duljine 5");

    uint16_t code = 0;
    for (char c : kmer)
    {
        code <<= 2; // pomakni ulijevo za 2 bita

        switch (c)
        {
        case 'A':
            code |= 0b00;
            break;
        case 'C':
            code |= 0b01;
            break;
        case 'G':
            code |= 0b10;
            break;
        case 'T':
            code |= 0b11;
            break;
        default:
            throw invalid_argument("Nevažeći znak u sekvenci");
        }
    }
    return code;
}

string decode(uint16_t code)
{
    string kmer(5, 'A'); // rezerviramo niz od 5 znakova, svi su 'A'

    for (int i = 4; i >= 0; --i)
    {
        uint16_t bits = code & 0b11; // uzmi zadnja 2 bita (najmanje značajne)
        switch (bits)
        {
        case 0b00:
            kmer[i] = 'A';
            break;
        case 0b01:
            kmer[i] = 'C';
            break;
        case 0b10:
            kmer[i] = 'G';
            break;
        case 0b11:
            kmer[i] = 'T';
            break;
        }
        code >>= 2; // pomakni udesno za 2 bita da dohvatiš sljedeće slovo
    }
    return kmer;
}

double euclid(const unordered_map<uint16_t, double> &dict)
{
    double sum = 0.0;

    for (const auto &pair : dict) {
        double freq = pair.second;
        sum += freq * freq;
    }

    return sqrt(sum);
}

int main()
{
    string path = "../Data/reading.fasta";

    // Kreiraj parser koji čita FASTA i sprema u objekte tipa Sequence
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

    // Parsiraj sve sekvence iz datoteke
    auto sekvence = parser->Parse(-1); // -1 znači: učitaj sve

    // auto &s = sekvence[0];
    // s->Print();

    // Testiranje enkodiranja i dekodiranja
    // string kmer = "ACGTG";
    // cout << "Originalni kmer: " << kmer << endl;
    // uint16_t encoded = encode(kmer);
    // cout << "Enkodirani kmer: " << encoded << endl;
    // string decoded = decode(encoded);
    // cout << "Dekodirani kmer: " << decoded << endl;

    for (const auto &s : sekvence)
    {
        string sekv = s->data;

        unordered_map<uint16_t, int> dict; // dict with kmer as key and its count in sequence as value

        for (int i = 0; i <= sekv.size() - 5; ++i)
        {
            string kmer = sekv.substr(i, 5);
            uint16_t encoded = encode(kmer);
            dict[encoded]++;
            // cout << "Kmer: " << kmer << endl;
        }

        int total_kmers = 0;
        for (const auto &pair : dict)
        {
            total_kmers += pair.second;
        }

        unordered_map<uint16_t, double> freq_dict; // dict with kmer as key and its relative count in sequence as value
        for (const auto &pair : dict)
        {
            double freq = static_cast<double>(pair.second) / total_kmers;
            freq_dict[pair.first] = freq;
            cout << "Kmer: " << decode(pair.first) << ", Relative count: " << freq << endl;
        }

        cout << "Euklidska norma vektora: " << euclid(freq_dict) << endl;
        cout << endl;
    }

    return 0;
}