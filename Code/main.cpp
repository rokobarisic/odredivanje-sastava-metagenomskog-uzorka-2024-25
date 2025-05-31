#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include "bioparser/fasta_parser.hpp"
#include <math.h>

#define KMER_LENGTH 5

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
    if (kmer.size() != KMER_LENGTH)
        throw invalid_argument("Sekvenca mora biti duljine " + to_string(KMER_LENGTH) + ", a tu je duljine " + to_string(kmer.size()) + kmer);

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
            throw invalid_argument("Nevažeći znak u sekvenci: " + string(1, c));
        }
    }
    return code;
}

string decode(uint16_t code)
{
    string kmer(KMER_LENGTH, 'A'); // rezerviramo niz od KMER_LENGTH znakova, svi su 'A'

    for (int i = KMER_LENGTH-1; i >= 0; --i)
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

    for (const auto &pair : dict)
    {
        double freq = pair.second;
        sum += freq * freq;
    }

    return sqrt(sum);
}

double scalar_product(const unordered_map<uint16_t, double> &dict1, const unordered_map<uint16_t, double> &dict2)
{
    double sum = 0.0;

    for (const auto &pair : dict1)
    {
        auto it = dict2.find(pair.first);
        if (it != dict2.end())
        {
            sum += pair.second * it->second;
        }
    }

    return sum;
}

unordered_map<uint16_t, double> get_freq_dict(const string &sekv)
{
    unordered_map<uint16_t, int> dict; // dict with kmer as key and its count in sequence as value

    for (int i = 0; i <= sekv.size() - KMER_LENGTH; i++)
    {
        if (sekv.size() < i + KMER_LENGTH)
            break;
        string kmer = sekv.substr(i, KMER_LENGTH);
        // if (kmer.size() < KMER_LENGTH) continue;
        if (kmer.find('N') != string::npos)
            continue;
        if (kmer.find('R') != string::npos)
            continue;
        if (kmer.find('Y') != string::npos)
            continue;
        if (kmer.find('K') != string::npos)
            continue;
        if (kmer.find('M') != string::npos)
            continue;
        if (kmer.find('S') != string::npos)
            continue;
        if (kmer.find('W') != string::npos)
            continue;
        if (kmer.find('B') != string::npos)
            continue;
        if (kmer.find('D') != string::npos)
            continue;
        if (kmer.find('H') != string::npos)
            continue;
        if (kmer.find('V') != string::npos)
            continue;
        if (kmer.find('X') != string::npos)
            continue;
        if (kmer.find('-') != string::npos)
            continue;
        uint16_t encoded = encode(kmer);
        dict[encoded]++;
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
    }

    return freq_dict;
}

int main()
{
    unordered_map<string, unordered_map<uint16_t, double>> reference_freq_dicts;
    unordered_map<string, int> readings_per_reference;

    // parsiramo referentne datoteke
    string path = "../Data/References";

    try
    {
        for (const auto &entry : fs::directory_iterator(path))
        {
            if (entry.path().filename() == ".gitkeep")
            {
                continue;
            }

            auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(entry.path().string());
            // cout << "Parsiram referentnu datoteku: " << entry.path().filename() << endl;
            auto sekvence = parser->Parse(-1);

            for (const auto &s : sekvence)
            {
                //readings_per_reference[s->id] = 0;
                string sekv = s->data;
                unordered_map<uint16_t, double> freq_dict = get_freq_dict(sekv);
                reference_freq_dicts[s->id] = freq_dict;
            }
        }

        // for (const auto &entry : reference_freq_dicts)
        // {
        //     cout << "Referentna datoteka: " << entry.first << endl;
        //     cout << "Euklidska norma vektora: " << euclid(entry.second) << endl;
        // }
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Greška: " << e.what() << endl;
    }

    string reading_path = "../Data/reading.fasta";

    // Kreiraj parser koji čita FASTA i sprema u objekte tipa Sequence
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(reading_path);

    // Parsiraj sve sekvence iz datoteke
    auto sekvence = parser->Parse(-1); // -1 znači: učitaj sve

    for (const auto &s : sekvence)
    {
        string sekv = s->data;

        unordered_map<uint16_t, double> freq_dict = get_freq_dict(sekv);

        // cout << "Euklidska norma vektora: " << euclid(freq_dict) << endl;
        // cout << endl;

        double max_similarity = 0.0;
        string most_similar_ref;

        for (const auto &entry : reference_freq_dicts)
        {
            const string &ref_name = entry.first;
            const unordered_map<uint16_t, double> &ref_freq_dict = entry.second;

            double reading_euclid = euclid(freq_dict);
            double reference_euclid = euclid(ref_freq_dict);
            double scalar_prod = scalar_product(freq_dict, ref_freq_dict);

            double similarity = scalar_prod / (reading_euclid * reference_euclid);
            //cout << "Slicnost sa referentnom datotekom " << ref_name << ": " << similarity << endl;
            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                most_similar_ref = ref_name;
            }
        }
        if (most_similar_ref.empty())
        {
            continue;
        }
        readings_per_reference[most_similar_ref]++;
        // cout << "ID: " << s->id << endl;
        // cout << "Najvise slici referentnom genomu: " << most_similar_ref << " (slicnost: " << max_similarity << ")" << endl << endl;
    }

    for (const auto &entry : readings_per_reference)
    {
        cout << "Referentna datoteka: " << entry.first << ", broj ocitanih sekvenci: " << entry.second << endl;
    }

    return 0;
}