#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <dirent.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>
using namespace std;


int main() {
    vector<vector <int> > vectorsq;
    cout << "build started" << endl;
    cout << "input: ";
    int terminput;
    cin >>  terminput;
    vector<int> input;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            input.push_back(i);
        }
    }
    for (int i = 0; i < vectorsq.size(); i++) {
        for (int j = 0; j < 4; j++) {
            cerr << vectorsq[i].at(j) << " ";
        }
        cerr << endl;
    }
    return 0;
}