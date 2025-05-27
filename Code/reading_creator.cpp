#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

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
            cout << "Datoteka ne postoji!" << endl;
            ofstream izlaz(path_reading);
            if (izlaz.is_open())
            {
                izlaz.close();
            }
        }
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Greška: " << e.what() << endl;
    }

    try
    {
        ofstream izlaz(path_reading);
        for (const auto &entry : fs::directory_iterator(path))
        {
            ifstream ocitanje(entry.path());
            int br_linija = 0;
            for (string linija; getline(ocitanje, linija);)
            {
                ++br_linija;
            }
            ocitanje.close();
            if (br_linija < 4)
            {
                continue;
            }
            int br_linija_kroz4 = br_linija / 4;

            int nr_readings = 1 + rand() % 10;

            for (int i = 0; i < nr_readings; i++)
            {
                ocitanje.open(entry.path());
                int random_number = rand() % (br_linija_kroz4 - 1);
                int index = random_number * 4 + 2;

                string linija;
                int trenutna = 1;
                while (getline(ocitanje, linija))
                {
                    if (trenutna == index - 1)
                    {
                        break;
                    }
                    else
                    {
                        trenutna++;
                    }
                }
                if (!linija.empty() && linija[0] == '@')
                {
                    linija[0] = '>';
                }
                izlaz << linija << endl;
                getline(ocitanje, linija);
                izlaz << linija << endl;
                ocitanje.close();
            }
        }
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Greška: " << e.what() << endl;
    }

    return 0;
}