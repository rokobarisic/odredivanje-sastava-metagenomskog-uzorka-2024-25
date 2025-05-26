#include <iostream>
#include <filesystem>
#include "bioparser/fasta_parser.hpp"

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

int main()
{
    string path = "../Data/reading.fasta"; // zamijeni sa stvarnim putem do FASTA datoteke

    // Kreiraj parser koji čita FASTA i sprema u objekte tipa Sequence
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

    // Parsiraj sve sekvence iz datoteke
    auto sekvence = parser->Parse(-1); // -1 znači: učitaj sve

    // Ispisi sve učitane sekvence
    for (const auto &s : sekvence)
    {
        s->Print();
    }

    return 0;
}