#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
double solveTime =0.;
double assembleTime =0.;
extern double calcstiffTime;
extern double contributeTime;

int main(int argc, char *argv[]) {
    calcstiffTime = 0.;
    contributeTime = 0.;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 0;
    pConfig.problem = "ESinSin";        //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Mixed";                    //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 4;                        //// How many refinements
    pConfig.debugger = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    
    if(pConfig.debugger)
    {
        std::string command = "cp ErroMixed.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
        FlushTable(pConfig,argv);
    }
    std::cout<<"calcstiffTime= "<< calcstiffTime << std::endl;

    return 0.;
}







