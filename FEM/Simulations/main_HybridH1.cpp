#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"
#include <TPZTimer.h>
#include "computStatist.h"
using namespace std;

#ifdef FEMCOMPARISON_TIMER
#include <math.h>//used only for timer statistics
#include <fstream>//used only for timer statistics
#include <vector>//used only for timer statistics

double solveTime =0.;
double assembleTime =0.;
extern double calcstiffTime;
extern double contributeTime; //  Total contribute time
double contributeTimeMaterial=0.;
double contributeTimeBoundary=0.;
double contributeTimeInterface=0.;
double interfaceTime=0.;
extern int64_t contributeCounter;
int64_t contributeMaterialCounter=0;
int64_t contributeBoundaryCounter=0;
double solveglobaltime;

vector<double> assembleTimeV, solveTimeV;
bool assembleTest=false;
bool solveTest=true;
int nThreads=0;
int nTestsAssemble=1;//number of tests for assemble
int nTestsSolve=5;//number of tests for assemble
//double mean(vector<double> vect);
//double standDeviation(vector<double>&);
//double CoefVariation(vector<double>&);
//void printTable(int dim,int nthreads, int refLevel, int test, bool MKL_FEMcomparison);
#endif

int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    TPZTimer timer;
    timer.start();
    calcstiffTime=0.;
    contributeTime =0;
    contributeCounter=0;
    bool atypical1=false;
    if((nTestsSolve>1) & (nTestsAssemble>1)){
        atypical1=true;
    }
    bool MKL_contribute;
#ifdef FEMCOMPARISON_USING_MKL
    MKL_contribute=true;
#endif
#endif
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    PreConfig pConfig;
    //Lectura por ficheros
       /*string line;
       ifstream myfile ("config.txt");
       getline (myfile,line);
       pConfig.k=stoi(line);
       getline (myfile,line);
       pConfig.n=stoi(line);
       getline (myfile,line);
       //pConfig.dim=stoi(line);
       //getline (myfile,line);
       pConfig.problem=line;
       getline (myfile,line);
       pConfig.approx=line;
       getline (myfile,line);
       pConfig.topology=line;
       getline (myfile,line);
       pConfig.refLevel=stoi(line);
       getline (myfile,line);
       myfile.close();*/
       //Fin Lectura Fichero

    pConfig.k = 1;//
    pConfig.n = 0;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Mixed";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 3;                       //// How many refinements
    pConfig.debugger = true;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel,pConfig,argv);
    DrawGeoMesh(config,pConfig);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.debugger){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
            FlushTable(pConfig,argv);
    }
    timer.stop();
    solveglobaltime = timer.seconds();
    if(assembleTest==true){
    cout<<"*********** Statistics for the assembly time ****************"<<endl;
    if(atypical1 == true)
        cout<<"Atypical experiment!!"<<endl;
    cout<<"reflevel: "<<pConfig.refLevel<<endl;
    if(MKL_contribute)
        cout<<"Using MKL in contributes: "<<"TRUE"<<endl;
    else
        cout<<"Using MKL in contributes: "<<"FALSE"<<endl;
    cout<<"Number of assembly threads: "<<nThreads<<endl;
    cout<<"Number of tests: "<<nTestsAssemble<<endl;
    cout<<"Average time(seconds): "<<mean(assembleTimeV)<<endl;
    cout<<"Coef. of variation: "<<100*CoefVariation(assembleTimeV)<<"%"<<endl;
    }
    if(solveTest==true){
    cout<<"*********** Statistics for the solve time ****************"<<endl;
    if(atypical1 == true)
        cout<<"Atypical experiment!!"<<endl;
    cout<<"reflevel: "<<pConfig.refLevel<<endl;
    if(MKL_contribute)
        cout<<"Using MKL in contributes: "<<"TRUE"<<endl;
    else
        cout<<"Using MKL in contributes: "<<"FALSE"<<endl;
    cout<<"Number of assembly threads: "<<nThreads<<endl;
    cout<<"Number of tests: "<<nTestsSolve<<endl;
    cout<<"Average time(seconds): "<<mean(solveTimeV)<<endl;
    cout<<"Coef. of variation: "<<100*CoefVariation(solveTimeV)<<"%"<<endl;
    }
/*
#ifdef FEMCOMPARISON_TIMER

     //Generar fichero
        string user = "ricardo-MAC";
        string experimentCode =pConfig.problem +"_"+pConfig.approx +"_"+pConfig.topology;
        ofstream archivo;
        archivo.open ("experiment.txt");
        archivo<<user<<endl;
        archivo<<experimentCode<<endl;
        archivo<<solveglobaltime<<endl;
        archivo<<solveTime<<endl;
        archivo<<assembleTime<<endl;
        archivo<<calcstiffTime<<endl;
        archivo<<contributeTime<<endl;
        archivo<<contributeTimeMaterial<<endl;
        archivo<<contributeTimeBoundary<<endl;
        archivo<<contributeTimeInterface<<endl;
        archivo<<pConfig.k<<endl;
        archivo<<pConfig.n<<endl;
        archivo<<pConfig.dim<<endl;
        archivo<<pConfig.refLevel<<endl;
#ifdef FEMCOMPARISON_USING_MKL
        archivo<<"MKL"<<endl;
#else
        archivo<<"NoMKL"<<endl;
#endif
        archivo.close();
    
    std::cout<<"solveglobaltime= "<< solveglobaltime << std::endl;
    std::cout<<"solveTime= "<< solveTime << std::endl;
    std::cout<<"assembleTime= "<< assembleTime << std::endl;
    std::cout<<"calcstiffTime= "<< calcstiffTime << std::endl;
    std::cout<<"contributeTimeInterface= "<< contributeTimeInterface << std::endl;
    std::cout<<"contributeTimeMaterial= "<< contributeTimeMaterial << std::endl;
    std::cout<<"contributeTimeBoundary= "<< contributeTimeBoundary << std::endl;
    //#ifdef OPTIMIZE_RUN_TIME
    std::cout<<"contributeTime= "<< contributeTime << std::endl;
    //#endif
    //cout<<"contadorTimeMaterial= "<<contadorTimeMaterial<<endl;
    cout<<"contributeCounter= "<<contributeCounter<<endl;
    cout<<"contributeMaterialCounter= "<<contributeMaterialCounter<<endl;
    cout<<"contributeBoundaryCounter= "<<contributeBoundaryCounter<<endl;
    
#endif
*/
    return 0;
}

/*#ifdef FEMCOMPARISON_TIMER2
double mean(vector<double> vect){
    vector<double>::iterator it;
    double sum=0;
    for(it=vect.begin(); it!=vect.end(); it++)
       sum+=*it;
    return sum/vect.size();
}

double standDeviation(vector<double>& vect){
    double sum=0;
    vector<double>::iterator it;
    double media=mean(vect);
    for(it=vect.begin(); it!=vect.end(); it++)
        sum+=pow(*it-media,2.);
    return pow(sum/vect.size(),0.5);
}

double CoefVariation(vector<double>& vect){
    double desvStand=standDeviation(vect);
    double media=mean(vect);

    return desvStand/media;
}


void printTable(int dim,int nthreads, int refLevel, int test, bool MKL_FEMcomparison)
{
    std::ofstream ofs ("Salida.txt", std::ofstream::out);
    ofs<<"H1Hybrid - 2D"<<endl;
    ofs<<"MKL contribute"<<    MKL_FEMcomparison<<endl;
    ofs<<"Number of threads"<<nthreads<<endl;
    ofs<<"Refinement level"<<refLevel<<endl;
    ofs<<"Number of tests"<<test<<endl;
    ofs.close();
}

#endif
*/




