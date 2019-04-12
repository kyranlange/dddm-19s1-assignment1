//
// Created by kyran on 24/03/19.
//

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "CAGenerator.h"

using namespace std;

bool debug = false;

void remove_carriage_return(std::string& line) {
    if (*line.rbegin() == '\r')     {
        line.erase(line.length() - 1);
    }
}

/**
 * Read the affinity matrix from a file
 * @param AffinityMatrixFileName The name of the file to read the matrix from.
 * @return
 */
vector<vector<int>> readAffinityMatrix(const std::string &AffinityMatrixFileName) {
    /*
        45 0 41 0
        0 71 1 71
        41 1 38 1
        0 71 1 71
     */
    vector<vector<int>> affinity;
    ifstream affinity_file(AffinityMatrixFileName);

    if (!affinity_file) {
        cerr << "Unable to open affinity file";
        exit(1);
    }

    string line;

    // work out number of attributes.
    int no_attributes = 0;
    while (std::getline(affinity_file, line))
        ++no_attributes;

    // Reset file
    affinity_file.clear();                 // clear fail and eof bits
    affinity_file.seekg(0, std::ios::beg); // back to the start!

    affinity = vector<vector<int>>(no_attributes, vector<int>(no_attributes, 0));

    int row = 0;
    while (getline(affinity_file, line)) {
        remove_carriage_return(line);
        istringstream iss(line);
        string affinity_line;
        getline(iss, affinity_line);

        size_t current;
        size_t next = -1;
        int column = 0;
        do {
            current = next + 1;
            next = affinity_line.find_first_of(' ', current);
            if (next != string::npos) {
                affinity[row][column] = stoi(affinity_line.substr(current, next - current));
            }
            column++;
        } while (next != string::npos);
        row++;
    }

    return affinity;
}

int bond (int attribute1, int attribute2, int index, const vector<vector<int>> &affinity, vector<int> &order) {
    if (attribute1 < 0 || attribute2 < 0) {
        return 0;
    }
    if (attribute1 > index || attribute2 > index) {
        return 0;
    }

    int sum = 0;

    for (int z = 0; z < affinity[0].size(); z++) {
        sum += affinity[z][order[attribute1]] * affinity[z][order[attribute2]];
    }

    return sum;
}

int contribution (int attribute1, int attribute2, int attribute3, const vector<vector<int>> &affinity, vector<int> &order) {
    int index = attribute2;
    return 2*bond(attribute1, attribute2, index, affinity, order)
           + 2*bond(attribute2, attribute3, index, affinity, order)
           - 2*bond(attribute1, attribute3, index, affinity, order);
}


vector<vector<int>> calculateClusteredAffinity(vector<vector<int>> affinity) {
    int no_attributes = affinity[0].size();
    vector<vector<int>> ca = vector<vector<int>>(no_attributes, vector<int>(no_attributes, 0));

    vector<int> order; // Order of attributes
    for (int i = 0; i < no_attributes; i++) {
        order.push_back(i);
    }

    ca[0] = affinity[0];
    ca[1] = affinity[1];

    int index = 2;

    while (index < no_attributes) { //Choose the 'best' location for attribute affinity[index]
        int max = -1;
        int loc = -1;
        int cont = 0;

        for (int i = 0; i < index; i++) {
            cont = contribution(i-1, index, i, affinity, order);
            if (cont > max) {
                max = cont;
                loc = i;
            }
        }

        cont = contribution(index - 1, index, index + 1, affinity, order);
        if (cont > max) {
            loc = index;
        }

        // Reshuffle Matrix
        for (int j = index; j > loc; j--) {
            ca[j] = ca[j-1];
            order[j] = order[j-1];
        }
        ca[loc] = affinity[index];
        order[loc] = index;

        if (debug) {
            cout << "Order: ";
            for (int i = 0; i < order.size(); i++) {
                cout << order[i] << ' ';
            }
            cout << endl;
        }

        index++;
    }

    if (debug) {
        cout << "Order: ";
        for (int i = 0; i < order.size(); i++) {
            cout << order[i] << ' ';
        }
        cout << endl;

        for (int i = 0; i < ca.size(); i++) {
            for (int j = 0; j < ca[i].size(); j++) {
                cout << ca[i][j] << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }

    // Shuffle columns
    for (int i = 0; i < no_attributes; i++) {
        for (int j = 0; j < no_attributes; j++) {
            ca[i][j] = affinity[order[i]][order[j]];
        }
    }

    return ca;
}



/**
 * Run the Bond Energy Algorithm to calculate the clustered affinity matrix CA of an Affinity Matrix
 * @param argc The number of arguments supplied
 * @param argv The arguments
 * @return 0 if the program is successful.
 */
int main (int argc, char *argv[]) {

    // Parse command line arguments
    std::string affinityMatrixFile;

    if (argc != 2) { // The first argument is the name of the program
        if (debug) {
            affinityMatrixFile = "aa_t.txt";
        } else {
            return -1;
        }
    } else {
        affinityMatrixFile = argv[1];
    }

    vector<vector<int>> affinity = readAffinityMatrix(affinityMatrixFile);

    if (debug) { //print affinity
        for (int i = 0; i < affinity.size(); i++) {
            for (int j = 0; j < affinity[i].size(); j++) {
                cout << affinity[i][j] << ' ';
            }
            cout << endl;
        }
    }

    vector<vector<int>> ca = calculateClusteredAffinity(affinity);

    // Print Clustered Affinity Matrix
    for (int i = 0; i < ca.size(); i++) {
        for (int j = 0; j < ca[i].size(); j++) {
            cout << ca[i][j] << ' ';
        }
        cout << endl;
    }
}