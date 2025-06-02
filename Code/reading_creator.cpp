#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdlib>
#include <ctime>

#define MIN_READINGS_PER_FILE 10000
#define MAX_READINGS_PER_FILE 100000

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

        ofstream exit(path_reading);
        if (!exit.is_open())
        {
            cerr << "Cannot open exit file!" << endl;
            return 1;
        }

        for (const auto &entry : fs::directory_iterator(path))
        {
            vector<string> lines;
            ifstream reading(entry.path());
            for (string line; getline(reading, line);)
            {
                lines.push_back(line);
            }
            reading.close();

            if (lines.size() < 4)
                continue;

            int no_lines_div4 = lines.size() / 4;
            int nr_readings = MIN_READINGS_PER_FILE + rand() % (MAX_READINGS_PER_FILE - MIN_READINGS_PER_FILE + 1); // random no of readings TBA to reading.fasta for this file

            for (int i = 0; i < nr_readings; ++i)
            {
                int random_number = rand() % no_lines_div4; 
                int index = random_number * 4 + 1; // random index of 4-line block TBA to reading.fasta

                string header = lines[index - 1];
                string sequence = lines[index];

                if (!header.empty() && header[0] == '@')
                    header[0] = '>'; // fasta id starts with '>', not '@'

                exit << header << endl;
                exit << sequence << endl;
            }
        }

        exit.close();
    }
    catch (const fs::filesystem_error &e)
    {
        cout << "Error: " << e.what() << endl;
    }

    return 0;
}
