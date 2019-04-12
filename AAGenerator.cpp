#include <utility>

#include <utility>

//
// Created by kyran on 24/03/19.
//

#include "AAGenerator.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

bool debug = true;

class attribute {
    std::string label;
public:
    attribute(string label, string name) {
        this->label = std::move(label);
        this->name = std::move(name);
    }

    std::string name;
};

class query {
    std::string label;
public:
    query(string label, string query) {
        this->label = std::move(label);
        this->query_s = std::move(query);
    }

    std::string query_s;
};

void remove_carriage_return(std::string& line) {
    if (*line.rbegin() == '\r')     {
        line.erase(line.length() - 1);
    }
}

std::vector<attribute> readAttributes(const std::string &attribute_file_name) {
 /* Label Name
    A1 PNO
    A2 PNAME
    A3 BUDGET
    A4 LOC
  */

    vector<attribute> attributes;
    ifstream attribute_file(attribute_file_name);

    if (!attribute_file) {
        cerr << "Unable to open attribute file";
        exit(1);
    }

    string line;
    while (getline(attribute_file, line)) {
        remove_carriage_return(line);
        istringstream iss(line);
        string label; string name;
        getline(iss, label, ' ');
        getline(iss, name, ' ');
        attributes.emplace_back(label, name);
    }
    attributes.erase(attributes.begin());
    return attributes;
}

std::vector<query> readQueries(const std::string &queries_file_name) {
    vector<query> queries;
    ifstream queries_file(queries_file_name);

    if (!queries_file) {
        cerr << "Unable to open query file";
        exit(1);
    }

    string line;
    while (getline(queries_file, line)) {
        remove_carriage_return(line);
        istringstream iss(line);
        string label; string query_s;
        getline(iss, label, ' ');
        getline(iss, query_s);
        queries.emplace_back(label, query_s);
    }

    return queries;
}

vector<vector<int>> readAccessFrequency(const std::string &freq_file_name, unsigned long no_queries) {
    /*
          S1 S2 S3
        q1 15 20 10
        q2 5 0 0
        q3 25 25 25
        q4 5 0 0
     */
    vector<vector<int>> access;
    ifstream access_file(freq_file_name);

    if (!access_file) {
        cerr << "Unable to open access file";
        exit(1);
    }

    string line;
    int query = -1;
    int sites = 0;
    while (getline(access_file, line)) {
        remove_carriage_return(line);
        istringstream iss(line);
        string label; string acc_f_line;
        getline(iss, label, ' ');
        getline(iss, acc_f_line);

        size_t current;
        size_t next = -1;
        int site = 0;
        do {
            current = next + 1;
            next = acc_f_line.find_first_of(' ', current);
            if (query < 0) {
                sites++;
                if (next == string::npos) {
                    access = vector<vector<int>>(no_queries, vector<int>(sites - 1, 0));
                }
            } else {
                access[query][site] = stoi(acc_f_line.substr(current, next - current));
                site++;
            }
            if (next == string::npos) {
                query++;
            }
        }
        while (next != string::npos);
    }

    return access;
}

vector<vector<int>> populateUsage(vector<attribute> attributes, vector<query> queries) {
    vector<vector<int>> usage = vector<vector<int>>(queries.size(), vector<int>(attributes.size(), 0));

    for (int i = 0; i < queries.size(); i++) {
        for (int j = 0; j < attributes.size(); j++) {
            if (queries[i].query_s.find(attributes[j].name) != std::string::npos) {
                usage[i][j] = 1;
            } else {
                usage[i][j] = 0;
            }
        }

    }

    if (debug) {
        cout << "Usage:" << endl;
        for (int i = 0; i < usage.size(); i++) {
            cout << "q" << i+1 << ' ';
            for (int j = 0; j < usage[i].size(); j++) {
                cout << usage[i][j] << ' ';
            }
            cout << endl << endl;
        }
    }

    return usage;
}

vector<vector<int>> populateAffinityMatrix(vector<vector<int>> usage, vector<vector<int>> access, unsigned long no_attributes, unsigned long no_queries) {
    vector<vector<int>> affinity = vector<vector<int>>(no_attributes, vector<int>(no_attributes, 0));
    int execution = 1;

    for (int i = 0; i < no_attributes; i++) {
        for (int j = 0; j < no_attributes; j++) { // for every attribute pair
            affinity[i][j] = 0; // initialise affinity to 0

            for (int k = 0; k < no_queries; k++) { // for each of the queries
                if (usage[k][i] == 1 && usage[k][j] == 1) { // if both of the attributes are used

                    for (int x = 0; x < access[k].size(); x++) { // sum up the usage across the sites
                        affinity[i][j] += execution * access[k][x];
                    }
                }
            }
        }
    }

    return affinity;
}

/**
 * Calculate the attribute affinity (AA) matrix given the input of files containing attributes, queries and
 * access frequencies.
 * @param argc The number of arguments supplied
 * @param argv The arguments
 * @return 0 if the program is successful.
 */
int main (int argc, char *argv[]) {

    // Parse command line arguments
    std::string attributes_file;
    std::string queries_file;
    std::string access_frequencies_file;

    if (argc != 4) { // The first argument is the name of the program
        if (debug) {
            attributes_file = "att_1.txt";
            queries_file = "query_1.txt";
            access_frequencies_file = "acc_1.txt";
        } else {
            return -1;
        }
    } else {
        attributes_file = argv[1];
        queries_file = argv[2];
        access_frequencies_file = argv[3];
    }

    vector<attribute> attributes = readAttributes(attributes_file);
    vector<query> queries = readQueries(queries_file);
    vector<vector<int>> usage = populateUsage(attributes, queries);
    vector<vector<int>> access = readAccessFrequency(access_frequencies_file, queries.size());

    vector<vector<int>> affinity = populateAffinityMatrix(usage, access, attributes.size(), queries.size());

    for (int i = 0; i < affinity.size(); i++) {
        for (int j = 0; j < affinity[i].size(); j++) {
            cout << affinity[i][j] << ' ';
        }
        cout << endl;
    }
}