#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"

#ifdef FEMCOMPARISON_TIMER
    double solveTime =0.;
    double assembleTime =0.;
    extern double calcstiffTime;
    extern double contributeTime;
    double contributeTimeMaterial=0.;
    double contributeTimeInterface=0;
#endif


int main(int argc, char *argv[]) {

#ifdef FEMCOMPARISON_TIMER
    calcstiffTime=0.;
    contributeTime =0;
#endif
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Tetrahedral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 1;                        //// How many refinements
    pConfig.debugger = true;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    DrawGeoMesh(config,pConfig);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.debugger){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
            FlushTable(pConfig,argv);
    }

#ifdef FEMCOMPARISON_TIMER
    std::cout<<"contributeTimeInterface= "<< contributeTimeInterface << std::endl;
    std::cout<<"contributeTimeMaterial= "<< contributeTimeMaterial << std::endl;
    std::cout<<"contributeTime= "<< contributeTime << std::endl;
#endif

    return 0;
}







