#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 1;
    pConfig.problem = "ESinSin";        //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Mixed";                    //// {"H1","Hybrid", "Mixed"}
    pConfig.refLevel = 1;                        //// How many refinements
    pConfig.debugger = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    
    std::string command = "cp ErroMixed.txt " + pConfig.plotfile + "/Erro.txt";
    system(command.c_str());
    FlushTable(pConfig,argv);    
    return 0.;
}







