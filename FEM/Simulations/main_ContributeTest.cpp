//
// Created by victor on 22/03/2021.
//

#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 5;                        //// How many refinements
    pConfig.debugger = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    // The following data is for convergence rate computations (if applicable)
    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    pConfig.hLog = pConfig.h*2;

    // Problem configuration and geometric mesh definition
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
//    DrawGeoMesh(config,pConfig);

    // The FEM abstractions and computations are here
    Solve(config,pConfig);

    // Flushing a .xls (spreadsheet) table containing the errors and convergence rates (if applicable)
    if(pConfig.debugger){
        FlushTable(pConfig,argv);
    }
    return 0;
}
