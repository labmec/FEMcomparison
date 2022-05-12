#include "ScriptTools.h"
#include "InputTreatment.h"

void Generate(PreConfig &pConfig){
        
    std::ofstream fileStream;
    pConfig.speedUpOfstream = &fileStream;
    //pConfig.speedUpOfstream->open(pConfig.automatedFileName +"/config.sh",std::ofstream::app);
    
    if( remove("config.sh") != 0) perror( "WARNING: FAIL TO REMOVE CONFIG.SH" );
    else puts( "CONFIG.SH SUCCESSFULLY CLEARED");
    pConfig.speedUpOfstream->open("config.sh",std::ofstream::app);

    fileStream << "#!/bin/bash\n";
    
    OfstreamPath(pConfig);
    
    for (int approxMethod=0; approxMethod < 3; approxMethod++){
        for (int topologyType=1; topologyType < 4; topologyType += 2){
            std::stringstream sufix;
            switch (topologyType) {
                case 0:
                    pConfig.topology = "Triangular";
                    pConfig.refLevel = pConfig.ref2D;
                    pConfig.dim = 2;
                    break;
                case 1:
                    pConfig.topology = "Quadrilateral";
                    pConfig.refLevel = pConfig.ref2D;
                    pConfig.dim = 2;
                    break;
                case 2:
                    pConfig.topology = "Tetrahedral";
                    pConfig.refLevel = pConfig.ref3D;
                    pConfig.dim = 3;
                    break;
                case 3:
                    pConfig.topology = "Hexahedral";
                    pConfig.refLevel = pConfig.ref3D;
                    pConfig.dim = 3;
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
            if (pConfig.target.automated)
                GenerateThreadSpan(pConfig);
            else if (pConfig.target.timeEfficiency)
                GenerateRefSpan(pConfig);
            else{
                DebugStop();
            }
        }
            
    }
    pConfig.speeedUpPath->close();
    pConfig.speedUpOfstream->close();
        
}

void Summarize(PreConfig &pConfig){
        
        std::ifstream *PathsIfs = new std::ifstream;
        std::ifstream *DataIfs = new std::ifstream;
        
        pConfig.stat.csv.resize(3);
        pConfig.stat.txt.resize(3);

        for (int i = 0; i < 3; i++){
            pConfig.stat.csv[i] = new std::ofstream;
            pConfig.stat.txt[i] = new std::ofstream;
        }

        PathsIfs->open("OfstreamPath.txt");
            std::string pathLine;
            int counter = 0;
            while (std::getline(*PathsIfs,pathLine))
            {
                counter++;
                std::string path;
                
                std::istringstream iss(pathLine);
                iss >> path;
                std::string fileName;
                if (pConfig.target.automated){
                    fileName = path + "/speedUp.csv";
                } else if (pConfig.target.timeEfficiency){
                    fileName = path + "/timeEfficiency.csv";
                }else{
                    DebugStop();
                }
                
                DataIfs->open(fileName);
                std::string dataLine;
                
                std::vector<std::string> avgName;
                avgName.push_back(path + "/assem");
                avgName.push_back(path + "/solve");
                avgName.push_back(path + "/total");

                if (avgName.size() != 3) {DebugStop();}
                for (int i = 0; i < 3; i++){
                    pConfig.stat.csv[i]->open(avgName[i]+"Avg.csv");
                    *pConfig.stat.csv[i] << "nThreads,avg,cvar";
                    std::string suffix = "";
                    std::string fileName;

                    if (pConfig.target.automated){
                        suffix = ",speedUp";
                        fileName = "SpeedUp.txt";
                    }else if (pConfig.target.timeEfficiency){
                        fileName = "timeEfficiency.txt";
                    }
                    *pConfig.stat.csv[i] << suffix << std::endl;

                    pConfig.stat.txt[i]->open(avgName[i] + fileName);
                    *pConfig.stat.txt[i] << "coordinates{" << std::endl;
                }
                
                int lineCounter = 0;
                int nSpan = -1;
                int loopSize = -1;
                int iterCounter = -1;
                bool iterEnd = true;
                
                int counter2 = 0;
                while (std::getline(*DataIfs,dataLine)){
            
                counter2++;
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
                    
                    ReadSimFile(cellVec,pConfig);
                                        
                    bool newSpanValue = false;
                    
                    if (pConfig.target.automated){
                        newSpanValue = (nSpan != pConfig.rSimulation.nThreads);
                    }else if (pConfig.target.timeEfficiency){
                        newSpanValue = (nSpan != pConfig.rSimulation.nRef);
                    }

                    if(newSpanValue && iterEnd){
                        if (pConfig.target.automated){
                            nSpan = pConfig.rSimulation.nThreads;
                        }else if (pConfig.target.timeEfficiency){
                            nSpan = pConfig.rSimulation.nRef;
                        }
                        loopSize = pConfig.rSimulation.iterNum;
                        
                        if (loopSize < 1)
                            break;
                        
                        iterCounter = 1;
                        iterEnd = false;
                        pConfig.stat.timeVec.Resize(loopSize,3);
                        
                    } else {
                        if (loopSize < 1)
                            break;
                        if (loopSize != 1 && pConfig.rSimulation.iterNum != -iterCounter){
                          std::cout << "loopsize " << loopSize << " pConfig.rSimulation.iterNum " << pConfig.rSimulation.iterNum << "iterCounter " << iterCounter << "iterEnd" << iterEnd << std::endl;
                                break;
                        }
                        if (pConfig.rSimulation.iterNum == 1-loopSize){
                            iterEnd = true;
                        } else {
                            
                            iterEnd = false;
                        }
                        iterCounter++;
                    }
                    
                    int ipos = (iterCounter-1);
                    
                    pConfig.stat.timeVec(ipos,0) = pConfig.rSimulation.assembleTime;
                    pConfig.stat.timeVec(ipos,1) = pConfig.rSimulation.solveTime;
                    pConfig.stat.timeVec(ipos,2) = pConfig.rSimulation.totalTime;
                        
                    if (iterEnd){
                        for (int i=0; i<3; i++){
                            pConfig.stat.avg[i]=0;
                            pConfig.stat.spu[i]=0;
                            pConfig.stat.cvar[i]=0;
                        }

                        for (int iter = 0; iter < loopSize; iter++){
                            for(int i=0; i < 3; i++)
                                pConfig.stat.avg[i] += pConfig.stat.timeVec(iter,i)/loopSize;
                        }
                        pConfig.stat.serialTime.Resize(3,0);
                        if (pConfig.target.automated){
                            if (pConfig.rSimulation.nThreads == 0){
                                for(int i=0; i < 3; i++){
                                    pConfig.stat.serialTime[i] = pConfig.stat.avg[i];
                                }
                            }
                            for (int i = 0; i < 3; i++){
                                pConfig.stat.spu[i] = pConfig.stat.serialTime[i]/pConfig.stat.avg[i];
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
                        
                        if (pConfig.target.automated){
                            for (int i = 0; i < 3; i++){
                                *pConfig.stat.csv[i] << pConfig.rSimulation.nThreads << "," << pConfig.stat.avg[i] << "," << pConfig.stat.cvar[i] << "," << pConfig.stat.spu[i] << std::endl;
                                *pConfig.stat.txt[i] << "(" << pConfig.rSimulation.nThreads << "," << pConfig.stat.spu[i] << ")\n";
                            }
                        } else if (pConfig.target.timeEfficiency){
                            for (int i = 0; i < 3; i++){
                                *pConfig.stat.csv[i] << pConfig.rSimulation.nRef << "," << pConfig.stat.avg[i] << "," << pConfig.stat.cvar[i] << std::endl;
                                *pConfig.stat.txt[i] << "(" << pConfig.stat.avg[i] << ",\n";
                            }
                        }
                        
                        pConfig.stat.timeVec.Resize(0,0);
                        loopSize = -1;
                        
                    }
                    lineCounter++;
                }
                
                DataIfs->close();
                for (int i = 0; i < 3; i++){
                    pConfig.stat.csv[i]->close();
                    *pConfig.stat.txt[i] << "}" << std::endl;
                    pConfig.stat.txt[i]->close();
                }
            }
}

void GenerateThreadSpan(PreConfig &pConfig){
    for (int nThreads=0; nThreads < pConfig.tData.maxThreads+1; nThreads+=2){
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

void GenerateRefSpan(PreConfig &pConfig){
    switch (pConfig.dim) {
        case 2:
            pConfig.tData.maxRef = pConfig.ref2D;
            break;
        case 3:
            pConfig.tData.maxRef = pConfig.ref3D;
            break;
        default:
            DebugStop();
            break;
    }
    for (int nref=2; nref < pConfig.tData.maxRef+1; nref++){
        for (int iterNum = 0; iterNum < pConfig.stat.nLoops; iterNum++){
            int iterCode = -999;
            if(iterNum == 0){
                iterCode = pConfig.stat.nLoops;
            }else {
                iterCode = -iterNum;
            }
            
            ManageOfstream(pConfig,nref,iterCode);

        }
    }
}
