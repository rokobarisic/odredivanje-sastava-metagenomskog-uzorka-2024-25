#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdlib>
#include <ctime>

using namespace std;
namespace fs = std::filesystem;

int main()
{
    srand(time(0));

    string path = "../Data/Readings";
    string path_reading = "../Data/reading.fasta";

    try
    {
        if (!fs::exists(path_reading))
        {
            ofstream init(path_reading);
            init.close();
        }

        ofstream izlaz(path_reading);
        if (!izlaz.is_open())
        {
            cerr << "Ne mogu otvoriti izlaznu datoteku!" << endl;
            return 1;
        }

        for (const auto &entry : fs::directory_iterator(path))
        {
            vector<string> linije;
            ifstream ocitanje(entry.path());
            for (string linija; getline(ocitanje, linija);)
            {
                linije.push_back(linija);
            }
            ocitanje.close();

            if (linije.size() < 4)
                continue;

            int br_linija_kroz4 = linije.size() / 4;
            int nr_readings = 10000 + rand() % 99000;

            for (int i = 0; i < nr_readings; ++i)
            {
                int random_number = rand() % br_linija_kroz4;
                int index = random_number * 4 + 1;

                string header = linije[index - 1];
                string sequence = linije[index];

                if (!header.empty() && header[0] == '@')
                    header[0] = '>';

                izlaz << header << endl;
                izlaz << sequence << endl;
            }
        }

        izlaz.close();
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Greška: " << e.what() << endl;
    }

    return 0;
}
