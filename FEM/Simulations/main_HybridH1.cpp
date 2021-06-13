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

#ifdef FEMCOMPARISON_TIMER
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
#endif
using namespace std;

int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    TPZTimer timer;
    timer.start();
    calcstiffTime=0.;
    contributeTime =0;
    contributeCounter=0;
#endif
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    PreConfig pConfig;
    //Lectura por ficheros
       string line;
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
       myfile.close();
       //Fin Lectura Fichero

    /*pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 5;*/                       //// How many refinements
    pConfig.debugger = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    //DrawGeoMesh(config,pConfig);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.debugger){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
            FlushTable(pConfig,argv);
    }
    timer.stop();
    solveglobaltime = timer.seconds();

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
    return 0;
}







