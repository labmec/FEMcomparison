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
#endif

int main(int argc, char *argv[]) {

    
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
    
    int maxThreads = 8;
    int ref2D = 8;
    int ref3D = 4;
    
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
    
    for (int approxMethod=0; approxMethod < 3; approxMethod++){
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
                pConfig.tData.nThreads = nthreads;
                pConfig.h = 1./pConfig.exp;
                ProblemConfig config;
                Configure(config,pConfig.refLevel,pConfig,argv);

                Solve(config,pConfig);

                pConfig.hLog = pConfig.h;
            }
            pConfig.speedUpOfstream->close();
        }
    }
     return 0;
}





