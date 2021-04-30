#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"
double solveTime =0.;
double assembleTime =0.;
extern double calcstiffTime;
extern double contributeTime; //  Total contribute time
double contributeTimeMaterial=0.;
double contributeTimeInterface=0.;

using namespace std;

int main(int argc, char *argv[]) {
    calcstiffTime=0.;
    contributeTime =0.;
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    PreConfig pConfig;
    //Lectura por ficheros
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
       //Fin LEctura Fichero
    
    
    
    //pConfig.k = 1;//
    //pConfig.n = 2;
    //pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst"}
    // pConfig.approx = "Hybrid";                   //// {"H1","Hybrid", "Mixed"}
    //pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    //pConfig.refLevel = 0;                        //// How many refinements
    //pConfig.debugger = true;                    //// Print geometric and computational mesh
    
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    DrawGeoMesh(config,pConfig);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.debugger){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
            FlushTable(pConfig,argv);
    }

    //Generar fichero
    string user = "romulo-MAC";
        string experimentCode =pConfig.problem +"_"+pConfig.approx +"_"+pConfig.topology;
        ofstream archivo;
        archivo.open ("experiment.txt");
        archivo<<experimentCode<<endl;
        archivo<<user<<endl;
        archivo<<solveTime<<endl;
        archivo<<assembleTime<<endl;
        archivo<<calcstiffTime<<endl;
        archivo<<contributeTime<<endl;
        archivo<<contributeTimeMaterial<<endl;
        archivo<<contributeTimeInterface<<endl;
        archivo<<pConfig.k<<endl;
        archivo<<pConfig.n<<endl;
        archivo<<pConfig.dim<<endl;
        archivo<<pConfig.refLevel<<endl;
        archivo<<"MKL"<<endl;
        archivo.close();
    return 0;
}








