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

// g++ main.cpp -std=c++17 -I../external/bioparser/include -o main -lz

// needed by Bioparser
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
        cout << "ID: " << id << endl;
        cout << "Sequence: " << data << endl << endl;
    }
};


// encode a kmer of length KMER_LENGTH into a 16-bit integer (works when KMER_LENGTH <= 8)
// input: kmer - string representing the kmer
// output: 16-bit integer representing the kmer
uint16_t encode(const std::string &kmer)
{
    if (kmer.size() != KMER_LENGTH)
        throw invalid_argument("The sequence must be of length " + to_string(KMER_LENGTH) + ", but here it is of length " + to_string(kmer.size()) + " " + kmer);

    uint16_t code = 0;
    for (char c : kmer)
    {
        code <<= 2; // shift left by 2 bits to make space for the next nucleotide

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
            throw invalid_argument("Invalid character in sequence: " + string(1, c));
        }
    }
    return code;
}

// decode a 16-bit integer into a kmer of length KMER_LENGTH
// input: code - 16-bit integer representing the kmer
// output: string representing the kmer
string decode(uint16_t code)
{
    string kmer(KMER_LENGTH, 'A');

    for (int i = KMER_LENGTH-1; i >= 0; --i)
    {
        uint16_t bits = code & 0b11;
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
        code >>= 2;
    }
    return kmer;
}

// calculate euclidean norm of a vector
// input: dict - unordered_map with kmer as key and its relative count in sequence as value
// output: euclidean norm of the vector
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

// calculate scalar product of two vectors
// input: dict1 and dict2 - unordered_maps with kmer as key and its relative count in sequence as value
// output: scalar product of the two vectors
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

// get frequency dictionary of a sequence
// input: seq - string representing the sequence
// output: unordered_map with kmer as key and its relative count in sequence as value (distribution vector)
unordered_map<uint16_t, double> get_freq_dict(const string &seq)
{
    unordered_map<uint16_t, int> dict; // dict with kmer as key and its count in sequence as value

    for (int i = 0; i <= seq.size() - KMER_LENGTH; i++)
    {
        if (seq.size() < i + KMER_LENGTH)
            break;
        string kmer = seq.substr(i, KMER_LENGTH);
        if (kmer.find('N') != string::npos) // skipping kmers with these characters because they are not nucleotides, but informative characters
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

    string path = "../Data/References";

    // parsing references
    try
    {
        for (const auto &entry : fs::directory_iterator(path))
        {
            if (entry.path().filename() == ".gitkeep")
            {
                continue;
            }

            auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(entry.path().string());
            auto sequences = parser->Parse(-1); // -1 means parse all sequences

            for (const auto &s : sequences)
            {
                string seq = s->data;
                unordered_map<uint16_t, double> freq_dict = get_freq_dict(seq);
                reference_freq_dicts[s->id] = freq_dict;
            }
        }
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Error: " << e.what() << endl;
    }

    string reading_path = "../Data/reading.fasta";

    // parsing readings
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(reading_path);

    auto sequences = parser->Parse(-1); 

    for (const auto &s : sequences)
    {
        string seq = s->data;

        unordered_map<uint16_t, double> freq_dict = get_freq_dict(seq);

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
    }

    string output_path = "../Data/out.txt";
    ofstream output_file(output_path);
    for (const auto &entry : readings_per_reference)
    {
        output_file << "Reference file: " << entry.first << ", number of read sequences: " << entry.second << endl;
    }

    output_file.close();

    return 0;
}