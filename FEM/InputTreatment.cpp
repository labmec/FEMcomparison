//
// Created by victor on 02/09/2020.
//

#include "InputTreatment.h"
#include "DataStructure.h"
#include "MeshInit.h"
#include "Tools.h"
#include <algorithm>

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
        config.k = atoi(argv[4]);
        config.n = atoi(argv[5]);
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
    if( remove( "Erro.txt" ) != 0) perror( "Error deleting file" );
    else puts( "Error log successfully deleted" );

    pConfig.Erro.open("Erro.txt",std::ofstream::app);
    pConfig.Erro << "----------COMPUTED ERRORS----------\n";

    pConfig.Log = new TPZVec<REAL>(pConfig.numErrors, -1);
    pConfig.rate = new TPZVec<REAL>(pConfig.numErrors, -1);

    ProblemConfig config;
    Configure(config,0,pConfig,argv);
    
    InitializeSpeedUp(pConfig);
    
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
    else
        puts( "Error log successfully deleted" );

    pConfig.timer.open(timer_name, std::ofstream::app);
}

void EvaluateEntry(int argc, char *argv[],PreConfig &pConfig){
    if(argc != 1 && argc != 8){
        std::cout << "Invalid entry";
        DebugStop();
    }
    if(argc == 8){
        pConfig.argc = argc;
        for(int i = 4; i < 8 ; i++)
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
        
        pConfig.topology = argv[3];
        pConfig.refLevel = atoi(argv[6]);
        pConfig.tData.nThreads = atoi(argv[7]);
        pConfig.k = atoi(argv[4]);
        pConfig.n = atoi(argv[5]);
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
    
    std::string resultsFile="SpeedUpResults";
    std::string command = "mkdir -p " + resultsFile;
    system(command.c_str());
    resultsFile += '/' +time;
    pConfig.speedUpFilePath = resultsFile+'/';
    
    command = "mkdir -p " + resultsFile;
    system(command.c_str());
    
}

