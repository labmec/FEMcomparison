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

vector<double> assembleTimeVec, solveTimeVec;
bool assembleTest=false;
bool solveTest=true;
int nThreads=0;
int nTestsAssemble=1;//number of tests for assemble
int nTestsSolve=5;//number of tests for assemble
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
    
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 4;                       //// How many refinements
    pConfig.debugger = false;                    //// Print geometric and computational mesh

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
    cout<<"Average time(seconds): "<<mean(assembleTimeVec)<<endl;
    //cout<<"Coef. of variation: "<<100*CoefVariation(assembleTimeVec)<<"%"<<endl;
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
    //cout<<"Average time(seconds): "<<mean(solveTimeVec)<<endl;
    //cout<<"Coef. of variation: "<<100*CoefVariation(solveTimeVec)<<"%"<<endl;
    }

    return 0;
}




