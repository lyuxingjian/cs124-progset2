#include<iostream>
#include<vector>
#include<cassert>
#include<string>
#include <fstream>
#include "matrix.h"
using namespace std;


int main(int argc, char *argv[]){
    (void) argc;

    string dstr(argv[2]), fname(argv[3]);
    ifstream infile(fname.c_str());

    int D = stoi(dstr);
    Matrix<int> X(D, D), Y(D, D);
    for (auto i=0; i<2*D; i++)
        for (auto j=0; j<D; j++)
            if (infile.is_open())
                infile >> (i<D?X(i, j):Y(i-D, j));
    auto Z = X.strassen(Y, 64);
    for (auto i=0; i<D; i++)
        cout << Z(i, i) << endl;
    infile.close();
}