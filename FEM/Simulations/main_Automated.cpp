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
#include <iostream>
using namespace std;

#ifdef FEMCOMPARISON_TIMER
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
        
        int maxThreads = pConfig.tData.maxThreads;
        
        std::ofstream fileStream;
        pConfig.speedUpOfstream = &fileStream;
        //pConfig.speedUpOfstream->open(pConfig.automatedFileName +"/config.sh",std::ofstream::app);
        
        if( remove("config.sh") != 0) perror( "WARNING: FAIL TO REMOVE CONFIG.SH" );
        else puts( "CONFIG.SH SUCCESSFULLY CLEARED");
        pConfig.speedUpOfstream->open("config.sh",std::ofstream::app);

        fileStream << "#!/bin/bash\n";
        
        OfstreamPath(pConfig);
        
        for (int approxMethod=0; approxMethod < 3; approxMethod++){
            for (int topologyType=0; topologyType<4; topologyType++){
                std::stringstream sufix;
                switch (topologyType) {
                    case 0:
                        pConfig.topology = "Triangular";
                        pConfig.refLevel = pConfig.ref2D;
                        break;
                    case 1:
                        pConfig.topology = "Quadrilateral";
                        pConfig.refLevel = pConfig.ref2D;
                        break;
                    case 2:
                        pConfig.topology = "Tetrahedral";
                        pConfig.refLevel = pConfig.ref3D;
                        break;
                    case 3:
                        pConfig.topology = "Hexahedral";
                        pConfig.refLevel = pConfig.ref3D;
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
                
                ManageOfstream(pConfig);

                for (int nThreads=0; nThreads < maxThreads+1; nThreads+=2){
                    for (int iterNum = 0; iterNum < pConfig.stat.nLoops; iterNum++){
                        int iterCode = -999;
                        if(iterNum == 0){
                            iterCode = pConfig.stat.nLoops;
                        }else {
                            iterCode = -iterNum;
                        }
                        
                        ManageOfstream(pConfig, nThreads,iterCode);

                    }
                }
            }
        }
        pConfig.speeedUpPath->close();
        pConfig.speedUpOfstream->close();
        
        
    } else if (pConfig.argc == 2){
        
        std::ifstream *PathsIfs = new std::ifstream;
        std::ifstream *DataIfs = new std::ifstream;
        
        pConfig.stat.csv.resize(3);
        for (int i = 0; i < 3; i++){
            pConfig.stat.csv[i] = new std::ofstream;
        }

        pConfig.stat.txt = new std::ofstream;
        
        PathsIfs->open("ofstreamPath.txt");
            std::string pathLine;
            while (std::getline(*PathsIfs,pathLine))
            {
                std::string path;
                
                std::istringstream iss(pathLine);
                iss >> path;
                std::string fileName = path + "/speedUp.csv";
                DataIfs->open(fileName);
                std::string dataLine;
                
                std::vector<std::string> avgName;
                avgName.push_back(path + "/assemAvg.csv");
                avgName.push_back(path + "/solveAvg.csv");
                avgName.push_back(path + "/totalAvg.csv");

                if (avgName.size() != 3) {DebugStop();}
                for (int i = 0; i < 3; i++){
                    pConfig.stat.csv[i]->open(avgName[i]);
                    *pConfig.stat.csv[i] << "nThreads,avg,cvar,speedUp" << std::endl;
                }
                
                pConfig.stat.txt->open(path + "/assemSpeedUp");
                *pConfig.stat.txt << "{" << std::endl;
                
                int lineCounter = 0;
                int nThreads = -1;
                int loopSize = -1;
                int iterCounter = -1;
                bool iterEnd = true;
                
                while (std::getline(*DataIfs,dataLine)){
                    if (lineCounter == 0){
                        lineCounter++;
                        continue;
                    }
                    std::vector<std::string> cellVec;
                    
                    std::stringstream str(dataLine);
                    std::string cell;
                    while(getline(str, cell,',')){
                        cellVec.push_back(cell);
                    }
                    
                    ReadAutomated(cellVec,pConfig);
                    
                    if(nThreads != pConfig.rAutomated.nThreads && iterEnd){
                        loopSize = pConfig.rAutomated.iterNum;
                        
                        if (loopSize < 1)
                            DebugStop();
                        
                        iterCounter = 1;
                        iterEnd = false;
                        pConfig.stat.timeVec.Resize(loopSize,3);
                        
                    } else {
                        if (loopSize < 1)
                            DebugStop();
                        if (pConfig.rAutomated.iterNum != -iterCounter)
                                DebugStop();
                        if (pConfig.rAutomated.iterNum == 1-loopSize){
                            iterEnd = true;
                        } else {
                            
                            iterEnd = false;
                        }
                        iterCounter++;
                    }
                    
                    int ipos = (iterCounter-1);
                    
                    pConfig.stat.timeVec(ipos,0) = pConfig.rAutomated.assembleTime;
                    pConfig.stat.timeVec(ipos,1) = pConfig.rAutomated.solveTime;
                    pConfig.stat.timeVec(ipos,2) = pConfig.rAutomated.totalTime;
                        
                    if (iterEnd){
                        
                        pConfig.stat.avg.Resize(3,0);
                        pConfig.stat.spu.Resize(3,0);
                        pConfig.stat.cvar.Resize(3,0);

                        for (int iter = 0; iter < loopSize; iter++){
                            for(int i=0; i < 3; i++)
                                pConfig.stat.avg[i] += pConfig.stat.timeVec(iter,i)/loopSize;
                        }
                        pConfig.stat.serialTime.Resize(3,0);
                        if (pConfig.rAutomated.nThreads == 0){
                            for(int i=0; i < 3; i++){
                                pConfig.stat.serialTime[i] = pConfig.stat.avg[i];
                            }
                        }
                        
                        for (int iter = 0; iter < loopSize; iter++){
                            for(int i=0; i < 3; i++){
                                pConfig.stat.cvar[i] += pow(pConfig.stat.timeVec(iter,i)-pConfig.stat.avg[i],2.);
                            }
                        }
                        
                        for(int i=0; i < 3; i++){
                            pConfig.stat.cvar[i] = pow(pConfig.stat.cvar[i]/loopSize,0.5);
                        }
            
                        for (int i = 0; i < 3; i++){
                            pConfig.stat.spu[i] = pConfig.stat.serialTime[i]/pConfig.stat.avg[i];
                        }
                        
                        for (int i = 0; i < 3; i++){
                            *pConfig.stat.csv[i] << pConfig.rAutomated.nThreads << "," << pConfig.stat.avg[i] << "," << pConfig.stat.cvar[i]<< "," <<pConfig.stat.spu[i] << std::endl;
                        }

                        *pConfig.stat.txt << "(" << pConfig.rAutomated.nThreads << "," << pConfig.stat.spu[0] << ")\n";

                        pConfig.stat.timeVec.Resize(0,0);
                        loopSize = -1;
                    }
                    
                    std::cout << std::endl;
                    lineCounter++;
                }
                DataIfs->close();
                std::cout << std::endl;
                
                for (int i = 0; i < 3; i++){
                    pConfig.stat.csv[i]->close();
                }
                
                pConfig.stat.txt->close();
            }
    }else{
        InitializeAutomated(pConfig);
        pConfig.exp *= pow(2,pConfig.refLevel);
        pConfig.h = 1./pConfig.exp;
        ProblemConfig config;
        Configure(config,pConfig.refLevel,pConfig,argv);
        Solve(config,pConfig);
    }
    
    return 0;
    
}





