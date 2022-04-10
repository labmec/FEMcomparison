#include<iostream>
#include <vector>
#include <math.h>
#include <fstream>
using namespace std;

double mean(vector<double> vect);
double standDeviation(vector<double>&);
double CoefVariation(vector<double>&);
//printTable(int dim,int nthreads, int refLevel, bool MKL_FEMcomparison);

inline double mean(vector<double> vect){
    vector<double>::iterator it;
    double sum=0;
    for(it=vect.begin(); it!=vect.end(); it++)
       sum+=*it;
    return sum/vect.size();
}

inline double standDeviation(vector<double>& vect){
    double sum=0;
    vector<double>::iterator it;
    double media=mean(vect);
    for(it=vect.begin(); it!=vect.end(); it++)
        sum+=pow(*it-media,2.);
    return pow(sum/vect.size(),0.5);
}

inline double CoefVariation(vector<double>& vect){
    double desvStand=standDeviation(vect);
    double media=mean(vect);

    return desvStand/media;
}


inline void printTable(int dim,int nthreads, int refLevel, int test, bool MKL_FEMcomparison)
{
    std::ofstream ofs ("Salida.txt", std::ofstream::out);
    ofs<<"H1Hybrid - 2D"<<endl;
    ofs<<"MKL contribute"<<    MKL_FEMcomparison<<endl;
    ofs<<"Number of threads"<<nthreads<<endl;
    ofs<<"Refinement level"<<refLevel<<endl;
    ofs<<"Number of tests"<<test<<endl;
    ofs.close();
}
