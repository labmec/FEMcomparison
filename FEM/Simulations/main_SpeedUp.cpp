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
//#include <chrono>
#include <math.h>//used only for timer statistics
#include <fstream>//used only for timer statistics
#include <vector>//used only for timer statistics

double solveTime =0.; //Time spent solving the linear system
double assembleTime =0.;
extern double calcstiffTime;
extern double contributeTime; //  Total contribute time
//extern int64_t contributeCounter;
//int64_t contributeMaterialCounter=0;
//int64_t contributeBoundaryCounter=0;
double solveglobaltime;
long long contributeTimeVol = 0;
long long contributeTimeBC = 0;
long long contributeTimeInterface = 0;
//vector<unsigned long int> contributeTimeVolVec; //Total volumetric contribute time
vector<unsigned long long> assembleTimeVec;
vector<unsigned long long> solveTimeVec;
vector<unsigned long long> contributeTimeVec;
vector<unsigned long long> contributeTimeBCVec;

bool contributeTest=true;// To activate the time measure of the three contributes
bool assembleTest=true;
bool solveTest=true;
int nThreads=0;//number of threads for assemble
int nTestsAssemble=1;//number of tests for assemble
int nTestsSolve=1;//number of tests for solving the system of equations
#endif
int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    contributeTime =0.;
    bool atypical1=false;
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
    pConfig.problem = "ESinSin";              //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";       //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.debugger = false;                  //// Print geometric and computational mesh
    pConfig.shouldColor =false;
    pConfig.isTBB = false;
    
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
    
    for (int approxMethod=0; approxMethod < 3; approxMethod++){
        int ref2D = 2, ref3D = 1, maxThreads = 8;
        for (int topologyType=0; topologyType<4; topologyType++){
            std::stringstream sufix;
            switch (topologyType) {
                case 0:
                    pConfig.topology = "Triangular";
                    pConfig.refLevel = ref2D;
                    break;
                case 1:
                    pConfig.topology = "Quadrilateral";
                    pConfig.refLevel = ref2D;
                    break;
                case 2:
                    pConfig.topology = "Tetrahedral";
                    pConfig.refLevel = ref3D;
                    break;
                case 3:
                    pConfig.topology = "Hexahedral";
                    pConfig.refLevel = ref3D;
                    break;
                default:
                    DebugStop();
                    break;
            }
            
            sufix << pConfig.topology << "_" << pConfig.approx << "_";
            switch (approxMethod) {
                case 0:
                    pConfig.approx = "Mixed";
                    pConfig.n=0;
                    sufix << "n0";
                    break;
                case 1:
                    pConfig.approx = "Mixed";
                    pConfig.n=1;
                    sufix << "n1";
                    break;
                case 2:
                    pConfig.approx = "Hybrid";
                    if(topologyType < 2){
                        pConfig.n=2;
                        sufix << "n2";
                    } else{
                        pConfig.n=3;
                        sufix << "n3";
                    }
                    break;
                default:
                    DebugStop();
                    break;
            }
            sufix << "_ref" << pConfig.refLevel;
            std::ofstream fileStream;
            pConfig.speedUpOfstream = &fileStream;
            pConfig.speedUpOfstream->open(pConfig.speedUpFilePath+sufix.str()+".csv",std::ofstream::app);
            *pConfig.speedUpOfstream << "threadNum," << "Assemble," << "Solver," << "Total\n";
            
            for (int nthreads=0; nthreads < maxThreads+1; nthreads+=2){
                nThreads = nthreads;
                pConfig.h = 1./pConfig.exp;
                ProblemConfig config;
                Configure(config,pConfig.refLevel,pConfig,argv);

                Solve(config,pConfig);

                pConfig.hLog = pConfig.h;
            }
            pConfig.speedUpOfstream->close();
        }
    }
    
    cout << "******* HybridH1 *******"<< endl;
    if(atypical1 == true)
        cout<<"Atypical experiment!!"<<endl;
    cout<<"reflevel: "<<pConfig.refLevel<<endl;
    if(MKL_contribute)
        cout<<"Using MKL in contributes: "<<"TRUE"<<endl;
    else
        cout<<"Using MKL in contributes: "<<"FALSE"<<endl;
    cout<<"Number of assembly threads: "<<nThreads<<endl;
#ifdef TIMER_CONTRIBUTE
    cout<<"*********** Statistics for the time of the three contributes *****"<<endl;
    cout<<"Number of assemble tests: "<<nTestsAssemble<<endl;
    cout<<"Average time(seconds): "<<mean(contributeTimeVec)*1e-9<<endl;
    cout<<"Coef. of variation: "<<100*CoefVariation(contributeTimeVec)<<"%"<<endl;
#endif
    cout<<"*********** Statistics for the assembly time *****"<<endl;
    cout<<"Number of assemble tests: "<<nTestsAssemble<<endl;
    cout<<"Average time(seconds): "<<mean(assembleTimeVec)*1E-9<<endl;
    cout<<"Coef. of variation: "<<100*CoefVariation(assembleTimeVec)<<"%"<<endl;
    vector<unsigned long long>::iterator it;
    for(it=assembleTimeVec.begin(); it!=assembleTimeVec.end(); it++)
        cout<<(*it)*1E-9<<endl;
    //}
    //if(solveTest==true){
        cout<<"*********** Statistics for the linear system solve time *****"<<endl;
            cout<<"Number of assembly threads: "<<nThreads<<endl;
            cout<<"Number of solve tests: "<<nTestsSolve<<endl;
            cout<<"Average time(seconds): "<<mean(solveTimeVec)*1E-9<<endl;
            cout<<"Coef. of variation: "<<100*CoefVariation(solveTimeVec)<<"%"<<endl;
    for(it=solveTimeVec.begin(); it!=solveTimeVec.end(); it++)
        cout<<(*it)*1E-9<<endl;
    //}
    printTableAssemble(pConfig.dim,MKL_contribute,pConfig.refLevel,nThreads, nTestsAssemble,assembleTimeVec);
    
     return 0;
}





