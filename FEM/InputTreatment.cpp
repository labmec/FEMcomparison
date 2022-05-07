//
// Created by victor on 02/09/2020.
//

#include "InputTreatment.h"
#include "DataStructure.h"
#include "MeshInit.h"
#include "Tools.h"
#include <algorithm>
#include <thread>

void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]){
    ReadEntry(config, pConfig);
    config.ndivisions = ndiv;
    config.dimension = pConfig.dim;
    config.prefine = false;
#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1;
    config.exact->fDimension = config.dimension;
#endif
    
    bool isOriginCentered = 0; /// Wheater the domain = [0,1]x^n or [-1,1]^n
    if(pConfig.type == 2) isOriginCentered = 1;

    TPZGeoMesh *gmesh;
    TPZManVector<int, 4> bcids(4, -1);
    gmesh = CreateGeoMesh(1, bcids, config.dimension,isOriginCentered,pConfig.topologyMode);

    UniformRefinement(config.ndivisions, gmesh);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    if (pConfig.type  == 2) {
        config.materialids.insert(2);
        config.materialids.insert(3);
        config.bcmaterialids.insert(-5);
        config.bcmaterialids.insert(-6);
        config.bcmaterialids.insert(-8);
        config.bcmaterialids.insert(-9);

        SetMultiPermeMaterials(config.gmesh);
    }

    if(pConfig.argc != 1) {
        config.k = atoi(argv[5]);
        config.n = atoi(argv[6]);
    }

    if(pConfig.debugger == true && ndiv != 0){
        DrawGeoMesh(config,pConfig);
    }
}

void ReadEntry(ProblemConfig &config, PreConfig &preConfig){

#ifndef OPTMIZE_RUN_TIME
    config.exact = new TLaplaceExample1;
    switch(preConfig.type){
        case 0:
            config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
            break;
        case 1:
            config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
            break;
        case 2:
            config.exact.operator*().fExact = TLaplaceExample1::ESteklovNonConst;
            preConfig.h*=2;
            break;
        default:
            DebugStop();
            break;
    }
#endif
    config.k = preConfig.k;
    config.n = preConfig.n;
    config.problemname = preConfig.problem;
}


void InitializeOutstream(PreConfig &pConfig, char *argv[]){
    //Error buffer
#ifdef USING_SPEEDUP
    InitializeSpeedUp(pConfig);
#endif
    if (pConfig.makeScript) InitializeSpeedUp(pConfig);
    if (pConfig.argc == 5 || pConfig.argc == 2) return;

    if( remove( "Erro.txt" ) != 0) perror( "Error deleting file" );
    //else puts( "Error log successfully deleted" );

    pConfig.Erro.open("Erro.txt",std::ofstream::app);
    pConfig.Erro << "----------COMPUTED ERRORS----------\n";

    pConfig.Log = new TPZVec<REAL>(pConfig.numErrors, -1);
    pConfig.rate = new TPZVec<REAL>(pConfig.numErrors, -1);

    ProblemConfig config;
    Configure(config,0,pConfig,argv);
        
    std::stringstream out;
    switch (pConfig.topologyMode) {
        case 1:
            pConfig.topologyFileName = "2D-Tri";
            break;
        case 2:
            pConfig.topologyFileName = "2D-Qua";
            break;
        case 3:
            pConfig.topologyFileName = "3D-Tetra";
            break;
        case 4:
            pConfig.topologyFileName = "3D-Hex";
            break;
        case 5:
            pConfig.topologyFileName = "3D-Prism";
            break;
        default:
            DebugStop();
            break;
    }

    switch(pConfig.mode) {
        case 0: //H1
            out << "H1_" <<  pConfig.topologyFileName << "_" << config.problemname << "_k-"
                << config.k;
            pConfig.plotfile = out.str();
            break;
        case 1: //Hybrid
            out << "Hybrid_" <<   pConfig.topologyFileName << "_" << config.problemname  << "_k-"
                << config.k << "_n-" << config.n;
            pConfig.plotfile = out.str();
            break;
        case 2: // Mixed
            out << "Mixed_" <<  pConfig.topologyFileName << "_" << config.problemname << "_k-"
                << config.k << "_n-" << config.n;
            pConfig.plotfile = out.str();
            break;
        default:
            std::cout << "Invalid mode number";
            DebugStop();
            break;
    }
    std::string command = "mkdir -p " + pConfig.plotfile;
    system(command.c_str());

    std::string timer_name = pConfig.plotfile + "/timer.txt";

    if( remove(timer_name.c_str()) != 0)
        perror( "Error deleting file" );
//    else
//        puts( "Error log successfully deleted" );

    pConfig.timer.open(timer_name, std::ofstream::app);
}

void EvaluateEntry(int argc, char *argv[],PreConfig &pConfig){
    pConfig.argc = argc;
    
    if (argc == 5){
        if (std::strcmp(argv[1],"generate") == 0){
            pConfig.makeScript = true;
            int maxnThreads = std::thread::hardware_concurrency();
            
            if (maxnThreads > 32) {
                std::cout << "maxnThreads = " << maxnThreads << "setting maxnThreads = 32\n";
                maxnThreads = 32;
            }
            pConfig.tData.maxThreads = maxnThreads;
            for(int i = 2; i < 5 ; i++)
                IsInteger(argv[i]);
            pConfig.ref2D = atoi(argv[2]);
            pConfig.ref3D = atoi(argv[3]);
            pConfig.stat.nLoops = atoi(argv[4]);
            
            if (pConfig.stat.nLoops < 1) {
                DebugStop();
            }
            return;
        }else {
            DebugStop();
        }
    }
    if (argc == 2){
        if (std::strcmp(argv[1],"summarize") == 0){
            return;
        } else {
            DebugStop();
        }
    }
    if(argc != 1 && argc != 11){
        std::cout << "Invalid entry";
        DebugStop();
    }
    if(argc == 11){
        for(int i = 5; i < 11 ; i++)
            IsInteger(argv[i]);
        if(std::strcmp(argv[2], "H1") == 0)
            pConfig.mode = 0;
        else if(std::strcmp(argv[2], "Hybrid") == 0) {
            pConfig.mode = 1;
        }
        else if(std::strcmp(argv[2], "Mixed") == 0) {
            pConfig.mode = 2;
            if(pConfig.n < 0) DebugStop();
        }
        else DebugStop();

        if(std::strcmp(argv[1], "ESinSin") == 0) {
            pConfig.type = 0;
            pConfig.problem = "ESinSin";
        }
        else if(std::strcmp(argv[1], "EArcTan") == 0) {
            pConfig.type = 1;
            if(pConfig.n < 1 ){
                std::cout << "Unstable method\n";
                DebugStop();
            }
            pConfig.problem = "EArcTan";
        }
        else if(std::strcmp(argv[1], "ESteklovNonConst") == 0) {
            pConfig.type = 2;
            if(pConfig.n < 0) DebugStop();
            pConfig.problem = "ESteklovNonConst";
        }
        else DebugStop();
        
        pConfig.approx=argv[2];
        pConfig.topology = argv[3];
        pConfig.automatedFileName = argv[4];
        pConfig.k = atoi(argv[5]);
        pConfig.n = atoi(argv[6]);
        pConfig.refLevel = atoi(argv[7]);
        pConfig.tData.nThreads = atoi(argv[8]);
        pConfig.tData.maxThreads = atoi(argv[9]);
        pConfig.stat.iterNum = atoi(argv[10]);

    }
    else{
        if (pConfig.approx == "H1") pConfig.mode = 0;
        else if (pConfig.approx == "Hybrid")  pConfig.mode = 1;
        else if (pConfig.approx == "Mixed") pConfig.mode = 2;
        else DebugStop();

        if (pConfig.problem== "ESinSin") pConfig.type= 0;
        else if (pConfig.problem=="EArcTan")  pConfig.type = 1;
        else if (pConfig.problem == "ESteklovNonConst") pConfig.type = 2;
        else DebugStop();
    }

    if (pConfig.topologyMode != -1) DebugStop();
    if (pConfig.topology == "Triangular") pConfig.topologyMode = 1;
    else if (pConfig.topology == "Quadrilateral") pConfig.topologyMode = 2;
    else if (pConfig.topology == "Tetrahedral") pConfig.topologyMode = 3;
    else if (pConfig.topology == "Hexahedral") pConfig.topologyMode = 4;
    else if (pConfig.topology == "Prism") pConfig.topologyMode = 5;
    if (pConfig.topologyMode == -1) DebugStop();

    if(pConfig.topologyMode < 3) pConfig.dim = 2;
    else pConfig.dim = 3;
    
    if (pConfig.topologyMode < 3 && pConfig.mode ==1 && pConfig.n < 2){
        DebugStop();
    }
    
    if (pConfig.topologyMode > 2 && pConfig.mode ==1 && pConfig.n < 3){
        DebugStop();
    }
}

void IsInteger(char *argv){
    std::istringstream ss(argv);
    int x;
    if (!(ss >> x)) {
        std::cerr << "Invalid number: " << argv << '\n';
        DebugStop();
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv << '\n';
        DebugStop();
    }
}

void IsInteger(std::string str){
    int n = str.length();
    // declaring character array
    char char_array[n + 1];
    // copying the contents of the
    // string to char array
    strcpy(char_array, str.c_str());
    IsInteger(char_array);
    return;
}


void isFloat( std::string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    if (!(iss.eof() && !iss.fail())) {DebugStop();}
    return;
}

void CharReplace(std::string &str, char find, char replace ) {
  std::replace(str.begin(), str.end(), find, replace);
}

void InitializeSpeedUp(PreConfig &pConfig){
    
    std::stringstream currentTime;
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    currentTime << std::ctime(&end_time);
    
    std::string time = currentTime.str();
    time.pop_back();
    CharReplace(time, ' ', '_');
    CharReplace(time, ':', '-');

    std::string resultsFile;
#ifdef USING_SPEEDUP
    resultsFile="SpeedUpResults";
    pConfig.speedUpFilePath = resultsFile+'/';
    DirectoryCreation(resultsFile,time);
#endif
    
    if (pConfig.makeScript){
        resultsFile = "AutomatedResults";
        pConfig.automatedFileName = resultsFile+ '/'+ time;
        DirectoryCreation(resultsFile,time);
    }
}

void InitializeAutomated(PreConfig &pConfig){
    std::stringstream fileNamess;

    std::string fileName = pConfig.automatedFileName;
    std::vector<std::string> elems = split(fileName, '/');
    DirectoryCreation(elems[0], elems[1]);
    
    fileNamess<< pConfig.approx << "n" << pConfig.n << "-" << pConfig.topology;
    fileName = fileNamess.str();
    pConfig.automatedFilePath = pConfig.automatedFileName + "/" + fileName;
    std::string command = "mkdir -p " + pConfig.automatedFilePath;
    system(command.c_str());
    
    pConfig.speedUpOfstream = new std::ofstream;
    fileNamess.str(std::string());
    fileNamess << pConfig.automatedFilePath << "/speedUp.csv";

    if (pConfig.tData.nThreads == 0 && pConfig.stat.iterNum > 0){
        std::string stringnom = fileNamess.str();
        if( remove(stringnom.c_str()) != 0) perror( "File successfully created" );
    }
    pConfig.speedUpOfstream = new std::ofstream;
    pConfig.speedUpOfstream->open(fileNamess.str(),std::ofstream::app);

    if (pConfig.tData.nThreads == 0 && pConfig.stat.iterNum > 0){
        *pConfig.speedUpOfstream << "nThreads,assembleTime,solveTime,totalTime,maxNthreads,iterNum" << std::endl;
    }
}

void DirectoryCreation(std::string &resultsFile, std::string &time){
    std::string command = "mkdir -p " + resultsFile;
    system(command.c_str());
    resultsFile += '/' +time;
    command = "mkdir -p " + resultsFile;
    system(command.c_str());
}

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


void OfstreamPath(PreConfig &pConfig){
    
    std::string filePath = "OfstreamPath.txt";

    //char* ccx = new char[filePath.length() + 1];
    //std::copy(filePath.begin(), filePath.end(), ccx);
    
    if( remove(filePath.c_str()) != 0) perror( "WARNING: FAIL TO REMOVE OFSTREAMPATH.TXT" );
    else puts( "OFSTREAMPATH.TXT SUCCESSFULLY CLEARED");

    pConfig.speeedUpPath = new std::ofstream;
    pConfig.speeedUpPath->open(filePath,std::ios::app);
    
}

void ManageOfstream(PreConfig &pConfig, int nThreads, int iterNum){

    *pConfig.speedUpOfstream << "./Automated " << pConfig.problem << " " << pConfig.approx << " " << pConfig.topology << " " << pConfig.automatedFileName <<" " << pConfig.k << " " << pConfig.n << " " << pConfig.refLevel << " " << nThreads << " " << pConfig.tData.maxThreads << " " << iterNum << std::endl;
}

void ManageOfstream(PreConfig &pConfig){

    std::stringstream fileNamess;
    fileNamess<< pConfig.approx << "n" << pConfig.n << "-" << pConfig.topology;
    pConfig.automatedFilePath = pConfig.automatedFileName + '/' + fileNamess.str() + '/';

    *pConfig.speeedUpPath << pConfig.automatedFilePath << std::endl;
}

void ReadAutomated(std::vector<std::string> &cellVec, PreConfig &pConfig){
    if (cellVec.size() != 6) {DebugStop(); }

    IsInteger(cellVec[0]);
    IsInteger(cellVec[4]);
    IsInteger(cellVec[5]);
    
    isFloat(cellVec[1]);
    isFloat(cellVec[2]);
    isFloat(cellVec[3]);

    pConfig.rAutomated.nThreads = atoi(cellVec[0].c_str());
    pConfig.rAutomated.assembleTime = atof(cellVec[1].c_str());
    pConfig.rAutomated.solveTime = atof(cellVec[2].c_str());
    pConfig.rAutomated.totalTime = atof(cellVec[3].c_str());
    pConfig.rAutomated.maxThreads = atoi(cellVec[4].c_str());
    pConfig.rAutomated.iterNum = atoi(cellVec[5].c_str());
    
}
