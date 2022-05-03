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
    pConfig.targetAutomated = true;
    pConfig.debugger = false;       // Print geometric and computational mesh
    pConfig.shouldColor =false;
    pConfig.isTBB = false;
    
    if (argc == 1){
        std::cout << "This program isn't supposed to run this way \n";
        DebugStop();
    }
    
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
    
    if (pConfig.makeScript){
        
        int maxThreads = pConfig.tData.nThreads;
        int ref2D = 2;
        int ref3D = 1;
        
        std::ofstream fileStream;
        pConfig.speedUpOfstream = &fileStream;
        pConfig.speedUpOfstream->open(pConfig.automatedFileName +"/config.sh",std::ofstream::app);
        fileStream << "#!/bin/bash\n";
        
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
                
                for (int nThreads=0; nThreads < maxThreads+1; nThreads+=2){
                    *pConfig.speedUpOfstream << "./Automated " << pConfig.problem << " " << pConfig.approx << " " << pConfig.topology << " " << pConfig.automatedFileName <<" " << pConfig.k << " " << pConfig.n << " " << pConfig.refLevel << " " << nThreads << std::endl;
                    
                }
            }
        }
        pConfig.speedUpOfstream->close();
    } else{
        InitializeAutomated(pConfig);
        pConfig.exp *= pow(2,pConfig.refLevel);
        pConfig.h = 1./pConfig.exp;
        ProblemConfig config;
        Configure(config,pConfig.refLevel,pConfig,argv);
        Solve(config,pConfig);
    }
    
    return 0;
    
}





