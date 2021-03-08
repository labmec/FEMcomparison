#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
double contributeTime=0.;
double assembleTime=0.;
double calcstiffGlobalTime=0.;
double solveGlobalTime=0.;
string timeSelector;
int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 2;//
    pConfig.n = 2;
    pConfig.dim = 2;
    pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Mixed";                   //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 4;                        //// How many refinements
    
    //LEctura por ficheros
    string line;
    ifstream myfile ("config.txt");
    getline (myfile,line);
    pConfig.k=stoi(line);
    getline (myfile,line);
    pConfig.n=stoi(line);
    getline (myfile,line);
    pConfig.dim=stoi(line);
    getline (myfile,line);
    pConfig.problem=line;
    getline (myfile,line);
    pConfig.approx=line;
    getline (myfile,line);
    pConfig.topology=line;
    getline (myfile,line);
    pConfig.refLevel=stoi(line);
    getline (myfile,line);
    myfile.close();
    timeSelector = line;
    //Fin LEctura Fichero
    
    pConfig.debugger = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel,pConfig,argv);
    DrawGeoMesh(config,pConfig);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.debugger){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
            FlushTable(pConfig,argv);
    }
    string user = "romulo-MAC";
    string experimentCode =pConfig.problem +"_"+pConfig.approx +"_"+pConfig.topology;
    ofstream archivo;
    archivo.open ("experiment.txt");
    archivo<<experimentCode<<endl;
    archivo<<user<<endl;
    archivo<<contributeTime<<endl;
    archivo<<assembleTime<<endl;
    archivo<<calcstiffGlobalTime<<endl;
    archivo<<solveGlobalTime<<endl;
    archivo<<pConfig.k<<endl;
    archivo<<pConfig.n<<endl;
    archivo<<pConfig.dim<<endl;
    archivo<<pConfig.refLevel<<endl;
    archivo<<timeSelector<<endl;
    archivo.close();
    return 0;
}







