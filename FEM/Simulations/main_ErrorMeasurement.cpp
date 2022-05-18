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
#include "ScriptTools.h"

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
    pConfig.target.errorMeasurement = true;
    pConfig.debugger = true;       // Print geometric and computational mesh
    pConfig.shouldColor =false;
    pConfig.isTBB = false;
    
    if (argc == 1){
        std::cout << "This program isn't supposed to run this way \n";
        DebugStop();
    }
    
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
    
    if (pConfig.makeScript){
        
        Generate(pConfig);
        
    }else if (pConfig.argc == 2){
        if (std::strcmp(argv[1],"simplify") == 0){
            std::vector<std::vector<int>*> dataPerSetup;
            dataPerSetup.resize(2);
            
            for (int i=0; i<2; i++){
                dataPerSetup[i] = new std::vector<int>();
            }

            CountDataPerSetup(dataPerSetup, pConfig);
            
            std::vector<std::ifstream*> PathsIfs;
            std::vector<std::ifstream*> DataIfs;
            PathsIfs.resize(2);
            
            int nErrors = 2;
            int nTimes = 3;
            std::vector<int> nMeasurements;
            nMeasurements.push_back(nErrors);
            nMeasurements.push_back(nTimes);
            DataIfs.resize(nErrors+nTimes);
            
            for (int i=0; i<2; i++){
                PathsIfs[i] = new std::ifstream;
            }
            for (int i=0; i<nErrors+nTimes; i++){
                DataIfs[i] = new std::ifstream;
            }
            
            std::vector<std::string> ofsPathFile;
            ofsPathFile.push_back("OfstreamPathError.txt");
            ofsPathFile.push_back("OfstreamPathTime.txt");
            
            for (int i=0; i<2; i++){
                PathsIfs[i]->open(ofsPathFile[i]);
            }
            
            //std::vector<std::vector<std::string>> containerName;
            std::vector<std::string> errorName, timerName, fileName;
            errorName.push_back("L2Error");
            errorName.push_back("EnergyError");
            timerName.push_back("assemtimeEfficiency");
            timerName.push_back("solvetimeEfficiency");
            timerName.push_back("totaltimeEfficiency");
//            containerName[0] = errorName;
//            containerName[1] = timerName;

            for (int i=0; i < errorName.size() ;i++){
                for (int j=0; j < timerName.size() ;j++){
                    fileName.push_back(errorName[i]+timerName[j]);
                }
            }
            
            for (int i=0; i < dataPerSetup[0]->size(); i++){
                std::vector<std::string> pathLine;
                pathLine.resize(2);
                for (int j=0; j < 2; j++){
                    std::getline(*(PathsIfs[j]),pathLine[j]);
                }
                for (int j=0; j< nMeasurements[0]; j++){
                    DataIfs[j]->open(pathLine[0]+errorName[j]+".txt");
                    
                    if(!DataIfs[j])
                       {
                           cout << "Error opening files!" << endl;
                           return 1;
                       }
                }
                
                for (int j=0; j< nMeasurements[1]; j++){
                    DataIfs[j+2]->open(pathLine[1]+timerName[j]+".txt");
                    if(!DataIfs[j])
                       {
                           cout << "Error opening files!" << endl;
                           return 1;
                       }
                }
                
                std::vector<std::ofstream*> DataOfs;
                DataOfs.resize(nMeasurements[0]*nMeasurements[1]);
                
                for (int i=0; i<nMeasurements[0]*nMeasurements[1]; i++){
                    DataOfs[i] = new std::ofstream;
                }
                
                for (int j=0; j< nMeasurements[0]; j++){
                    for (int k=0; k < nMeasurements[1]; k++){
                        DataOfs[3*j+k]->open(pathLine[0]+errorName[j]+timerName[k]+".txt");
                    }
                }

                for (int j=0; j < (*(dataPerSetup[0]))[i]; j++){
                    std::vector<std::string> rTime, rError, wVar;
                    rTime.resize(nMeasurements[1]);
                    rError.resize(nMeasurements[0]);
                    wVar.resize(nMeasurements[0]*nMeasurements[1]);

                    for (int k=0; k < 3; k++){
                        *DataIfs[k+2] >> rTime[k];
                        if (!(j == 0 || j == (*(dataPerSetup[0]))[i]-1)){
                            rTime[k].pop_back();
                        }

                    }
                    
                    for (int k=0; k < 2; k++){
                        *DataIfs[k] >> rError[k];
                        if (j == 0 || j == (*(dataPerSetup[0]))[i]-1){
                            rError[k] = "";
                        }
                    }
                    
                    for (int k=0; k < 2; k++){
                        for (int z=0; z < 3; z++){
                            wVar[3*k+z] = rTime[z]+rError[k];
                            *DataOfs[3*k+z] << wVar[3*k+z] << "\n";
                        }
                    }
                }

                
                for (int j=0; j< nMeasurements[0]; j++){
                    for (int k=0; k < nMeasurements[1]; k++){
                        DataOfs[3*j+k]->close();
                    }
                }
                for (int j=0; j<  nMeasurements[0]+nMeasurements[1]; j++){
                    (DataIfs[j])->close();
                }
            }
            
            for (int i=0; i<2; i++){
                PathsIfs[i]->close();
            }
            
            
        }else
            
            Summarize(pConfig);
        
    }else{
        InitializeOfsSim(pConfig);
        pConfig.exp *= pow(2,pConfig.refLevel);
        pConfig.h = 1./pConfig.exp;
        ProblemConfig config;
        Configure(config,pConfig.refLevel,pConfig,argv);
        Solve(config,pConfig);
    }
    
    return 0;
    
}





