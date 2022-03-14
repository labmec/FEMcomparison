//
// Created by victor on 02/09/2020.
//

#include "Solver.h"
#include <TPZMultiphysicsCompMesh.h>
#include "pzanalysis.h"
#include "DataStructure.h"
#include "MeshInit.h"
#include "TPZCompMeshTools.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "Tools.h"
#include "Output.h"
#include "pzvisualmatrix.h"
#include "MeshInit.h"
#include "TPZTimer.h"
#include <chrono>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include "omp.h"

#ifdef PZ_LOG
static TPZLogger loggerST("solveTime");
static TPZLogger loggerAT("assembleTime");
#endif

#ifdef FEMCOMPARISON_TIMER
    extern std::vector<unsigned long long> assembleTimeVec;
    extern std::vector<unsigned long long> solveTimeVec;
    extern std::vector<unsigned long long> contributeTimeVec;
    extern std::vector<unsigned long long> contributeTimeBCVec;
    extern long long contributeTimeVol;
    extern long long contributeTimeBC;
    extern long long contributeTimeInterface;
#endif

void Solve(ProblemConfig &config, PreConfig &preConfig){

#ifdef FEMCOMPARISON_TIMER
    extern double solveTime;
    extern int nTestsSolve;
#endif
    TPZCompMesh *cmesh = InsertCMeshH1(config,preConfig);
    TPZMultiphysicsCompMesh *multiCmesh = new TPZMultiphysicsCompMesh(config.gmesh);
    int interfaceMatID = -10;
    int hybridLevel = 1;

    const clock_t start = clock();

    switch(preConfig.mode){
        case 0: //H1
            TPZCompMeshTools::CreatedCondensedElements(cmesh, false, false);
            SolveH1Problem(cmesh, config, preConfig);
            break;
        case 1: //Hybrid
            CreateHybridH1ComputationalMesh(multiCmesh, interfaceMatID,preConfig, config,hybridLevel);
            for(int ii=0;ii<1;ii++){
                //SolveHybridH1Problem(multiCmesh, interfaceMatID, config, preConfig,hybridLevel);
                int nthread = 6;
                //calcstiffTestOMP(multiCmesh,nthread);
                calcstiffTestTBB(multiCmesh,6);
            }
            break;
        case 2: //Mixed
            CreateMixedComputationalMesh(multiCmesh, preConfig, config);
            {
#ifdef PZ_LOG
                TPZTimer timer;
                if(loggerST.isDebugEnabled())
                timer.start();
#endif
            SolveMixedProblem(multiCmesh, config, preConfig);
#ifdef PZ_LOG
                if(loggerST.isDebugEnabled()){
                timer.stop();
                solveTime+=timer.seconds();
                }
#endif
            }
            break;
            default:
            DebugStop();
            break;
    }
    FlushTime(preConfig,start);

    if(preConfig.debugger) DrawCompMesh(config,preConfig,cmesh,multiCmesh);
}

void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh) {

    std::stringstream ref;
    ref << "_ref-" << 1/preConfig.h <<" x " << 1/preConfig.h;
    std::string refinement =  ref.str();

    std::ofstream out(preConfig.plotfile + "/gmesh"+ refinement + ".vtk");
    std::ofstream out2(preConfig.plotfile + "/gmesh"+ refinement + "txt");
    std::ofstream out3(preConfig.plotfile + "/cmesh.txt");

    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
    config.gmesh->Print(out2);

    if (preConfig.mode == 0) cmesh->Print(out3);
    else multiCmesh->Print(out3);
}

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, PreConfig &pConfig, ProblemConfig &config){
    InsertMaterialMixed(cmesh_Mixed, config,pConfig);
    TPZManVector<TPZCompMesh *, 4> meshvector(4);
    CreateMixedAtomicMeshes(meshvector,pConfig, config);

    TPZManVector<int> active(4, 1);

    cmesh_Mixed->BuildMultiphysicsSpace(active,meshvector);
    CreateCondensedMixedElements(cmesh_Mixed);
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
}
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed){

    int numActiveSpaces = cmesh_Mixed->MeshVector().size();
    if(numActiveSpaces == 4){
        int64_t nconnects = cmesh_Mixed->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = cmesh_Mixed->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == 4) c.IncrementElConnected();
        }
    }
    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh_Mixed, keeponelagrangian, keepmatrix);
}

void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int &interFaceMatID , PreConfig &pConfig, ProblemConfig &config,int hybridLevel){
    auto spaceType = TPZCreateMultiphysicsSpace::EH1Hybrid;
    if(hybridLevel == 2) {
        spaceType = TPZCreateMultiphysicsSpace::EH1HybridSquared;
    }
    else if(hybridLevel != 1) {
        DebugStop();
    }

    TPZCreateMultiphysicsSpace createspace(config.gmesh, spaceType);
    //TPZCreateMultiphysicsSpace createspace(config.gmesh);
    std::cout << cmesh_H1Hybrid->NEquations();

    createspace.SetMaterialIds({1,2,3}, {-6,-5,-2,-1});
    createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    int pOrder = config.n+config.k;
    createspace.CreateAtomicMeshes(meshvec,pOrder,config.k);

    InsertMaterialHybrid(cmesh_H1Hybrid, config,pConfig);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);

    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    interFaceMatID = createspace.fH1Hybrid.fLagrangeMatid.first;

}

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &pConfig){

#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = -1;
#endif
    
    std::cout << "Solving H1 " << std::endl;

    TPZAnalysis an(cmeshH1);

#ifdef FEMCOMPARISON_USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    std::set<int> matids;
    for (auto matid : config.materialids) matids.insert(matid);

    for(auto mat:config.bcmaterialids){
        matids.insert(mat);
    }

    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui

    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);

    ////Calculo do erro
    std::cout << "Computing Error H1 " << std::endl;

#ifndef OPTMIZE_RUN_TIME
    an.SetExact(config.exact.operator*().ExactSolution());
#endif
    
    StockErrorsH1(an,cmeshH1,pConfig.Erro,pConfig.Log,pConfig);

    ////PostProcess
    if(pConfig.debugger) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        vecnames.Push("Derivative");
        scalnames.Push("ExactSolution");

        int dim = cmeshH1->Reference()->Dimension();

        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution=0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution,dim);
    }
    std::cout << "FINISHED!" << std::endl;
}

void calcstiffTestOMP(TPZCompMesh *cmesh,int nthread){
    auto beginCalcStiff = std::chrono::high_resolution_clock::now();
    int64_t nelem = cmesh->NElements();
    omp_set_num_threads(nthread);

#pragma omp parallel for
    for (int64_t iel = 0; iel < nelem; iel++)
    {
        TPZCompEl *el = cmesh->Element(iel);
        if (!el) continue;
        TPZElementMatrix ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);

        el->CalcStiff(ek, ef);
    }
    auto endCalcStiff = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(endCalcStiff - beginCalcStiff);
    unsigned long int duration = static_cast<unsigned long int>(elapsed.count());
    std::cout << " CalcStiff duration omp = " << duration*1E-9 << " seconds" << std::endl;
}

void calcstiffTestTBB(TPZCompMesh *cmesh,int nthread){
    auto beginCalcStiff = std::chrono::high_resolution_clock::now();

    int64_t nelem = cmesh->NElements();
    tbb::task_scheduler_init init(nthread);
    tbb::parallel_for( tbb::blocked_range<int64_t>(0,nelem),
                      [&](tbb::blocked_range<int64_t> r){
    for (int64_t iel = r.begin(); iel < r.end(); iel++)
    {
        TPZCompEl *el = cmesh->Element(iel);
        if (!el) continue;
        TPZElementMatrix ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);

        el->CalcStiff(ek, ef);
    }
    });
    auto endCalcStiff = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(endCalcStiff - beginCalcStiff);
    unsigned long int duration = static_cast<unsigned long int>(elapsed.count());
    std::cout << " CalcStiff duration tbb= " << duration*1E-9 << " seconds" << std::endl;
}

                      
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId, struct ProblemConfig config,struct PreConfig &pConfig,int hybridLevel){
#ifdef FEMCOMPARISON_TIMER
    extern double solveTime;
    extern double assembleTime;
    extern int nTestsAssemble;
    extern int nTestsSolve;
#endif
#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1;
#endif
    
    std::cout << "Solving HYBRID_H1 " << std::endl;

    TPZAnalysis an(cmesh_H1Hybrid);

    extern int nThreads;

//    {
//            TPZCompMesh* fMesh = cmesh_H1Hybrid;
//            int64_t neq = fMesh->NEquations();
//
//            TPZFMatrix<STATE> *stiff = new TPZFMatrix<STATE>(neq,neq);//new TPZFYsmpMatrix(neq,neq);
//
//            TPZFMatrix<STATE> rhs;
//            rhs.Redim(neq,1);
//
//            int64_t iel;
//            int64_t nelem = fMesh->NElements();
//            TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//
//            TPZAdmChunkVector<TPZCompEl*> &elementvec = fMesh->ElementVec();
//
//            std::map<int ,TPZMaterial*> matvec = fMesh->MaterialVec();
//            std::set<int> fMaterialIds;
//            for(std::map<int ,TPZMaterial*>::iterator it = matvec.begin(); it != matvec.end(); ++it) {
//                fMaterialIds.insert(it->first);
//
//            }
//
//        int64_t count = 0;
//        auto beginCalcStiff = std::chrono::high_resolution_clock::now();
//
//        for (iel = 0; iel < nelem; iel++)
//        {
//                TPZCompEl *el = elementvec[iel];
//                if (!el) continue;
//                int matid = 0;
//                TPZGeoEl *gel = el->Reference();
//                if (gel) {
//                    matid = gel->MaterialId();
//                }
//                int matidsize = fMaterialIds.size();
//                if(matidsize){
//                    if(!el->NeedsComputing(fMaterialIds)) continue;
//                }
//
//                count++;
//                if (!(count % 1000)) {
//                    std::cout << '*';
//                    std::cout.flush();
//                }
//                if (!(count % 20000)) {
//                    std::cout << "\n";
//                }
//
//                //calcstiff.start();
//                ek.Reset();
//                ef.Reset();
//
//                el->CalcStiff(ek,ef);
//
//                //calcstiff.stop();
//                //assemble.start();
//                /*
//                if (!ek.HasDependency()) {
//                    ek.ComputeDestinationIndices();
//                    fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//
//                    stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//
//                    rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//                              test2.Print("matriz de rigidez interface",std::cout);
//
//                } else {
//                    // the element has dependent nodes
//                    ek.ApplyConstraints();
//                    ef.ApplyConstraints();
//                    ek.ComputeDestinationIndices();
//                    fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                    stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//                    rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//
//                }*/
//
//                //assemble.stop();
//            }
//        auto endCalcStiff = std::chrono::high_resolution_clock::now();
//        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(endCalcStiff - beginCalcStiff);
//        unsigned long int duration = static_cast<unsigned long int>(elapsed.count());
//        std::cout << " CalcStiff duration = " << duration*1E-9 << " seconds" << std::endl;
//
//    }
        //return 0;
    
#ifdef PZ_USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(nThreads);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    //    strmat.SetNumThreads(2);
    //    strmat.SetDecomposeType(ELDLt);
    
    
    TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
#endif
    std::set<int> matIds;
    for (auto matid : config.materialids) matIds.insert(matid);
    for (auto matidbc : config.bcmaterialids) matIds.insert(matidbc);

    matIds.insert(InterfaceMatId);
    strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
#ifdef FEMCOMPARISON_TIMER
    for(int i=0;i<nTestsAssemble;i++){
        auto beginAss = std::chrono::high_resolution_clock::now();
#ifdef TIMER_CONTRIBUTE
        contributeTimeVol = 0;
        contributeTimeBC = 0;
        contributeTimeInterface = 0;
#endif
#endif
        an.Assemble();
#ifdef FEMCOMPARISON_TIMER
#ifdef TIMER_CONTRIBUTE
        contributeTimeVec.push_back(contributeTimeVol+contributeTimeBC+contributeTimeInterface);
        //std::cout<<contributeTimeInterface<<std::endl;
        //contributeTimeBCVec.push_back(contributeTimeBC);
#endif
        auto endAss = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(endAss - beginAss);
        assembleTimeVec.push_back(static_cast<unsigned long int>(elapsed.count()));
        }
        //assembleTimeVec.push_back(static_cast<double>(timer.seconds()));
#endif
    //return 0;
#ifdef FEMCOMPARISON_TIMER
    for(int i=0;i<nTestsSolve;i++){
        auto begin = std::chrono::high_resolution_clock::now();
#endif
        an.Solve();
#ifdef FEMCOMPARISON_TIMER
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedSol = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        solveTimeVec.push_back(static_cast<unsigned long int>(elapsedSol.count()));
        //solveTimeVec.push_back(static_cast<double>(timer.seconds()));
    }
#endif

    int64_t nelem = cmesh_H1Hybrid->NElements();
    cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
    cmesh_H1Hybrid->ExpandSolution();
    cmesh_H1Hybrid->ElementSolution().Redim(nelem, 5);

    if(pConfig.debugger) {
        std::cout << "Computing Error HYBRID_H1 " << std::endl;
#ifndef OPTMIZE_RUN_TIME
        an.SetExact(config.exact.operator*().ExactSolution());
#endif
        ////Calculo do erro
        StockErrors(an,cmesh_H1Hybrid,pConfig.Erro,pConfig.Log,pConfig);
        std::cout << "DOF = " << cmesh_H1Hybrid->NEquations() << std::endl;
        ////PostProcess
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("PressureExact");
        vecnames.Push("Flux");

        int dim = pConfig.dim;
        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution = 3;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
    
}

using namespace std;

void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &pConfig) {
#ifdef FEMCOMPARISON_TIMER
    extern double solveTime;
    extern double assembleTime;
    extern int nTestsAssemble;
    extern int nTestsSolve;
#endif
#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1;
#endif
    bool optBW = true;

    std::cout << "Solving Mixed " << std::endl;
    TPZAnalysis an(cmesh_Mixed, optBW); //Cria objeto de análise que gerenciará a analise do problema
    if(false){
        cout<<"Total ecuaciones:"<<an.Solution().Rows()<<endl;
    }
    extern int nThreads;
    //MKL solver
#ifdef FEMCOMPARISON_USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_Mixed);
    //strmat.SetNumThreads(8);
    strmat.SetNumThreads(nThreads);
#else
    TPZSkylineStructMatrix strmat(cmesh_Mixed);
    strmat.SetNumThreads(0);
#endif
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    
#ifdef FEMCOMPARISON_TIMER
    for(int i=0;i<nTestsAssemble;i++){
        auto begin = std::chrono::high_resolution_clock::now();
#ifdef TIMER_CONTRIBUTE
        contributeTimeVol = 0;
        contributeTimeBC = 0;
        contributeTimeInterface = 0;
#endif
#endif
        an.Assemble();
#ifdef FEMCOMPARISON_TIMER
#ifdef TIMER_CONTRIBUTE
        contributeTimeVec.push_back(contributeTimeVol+contributeTimeBC+contributeTimeInterface);
        //std::cout<<contributeTimeInterface<<std::endl;
        //contributeTimeBCVec.push_back(contributeTimeBC);
#endif
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        assembleTimeVec.push_back(static_cast<unsigned long int>(elapsed.count()));
        //assembleTimeVec.push_back(static_cast<double>(timer.seconds()));
    }
#endif
    //return 0;
#ifdef FEMCOMPARISON_TIMER
    for(int i=0;i<nTestsSolve;i++){
        auto begin = std::chrono::high_resolution_clock::now();
#endif
        an.Solve();
#ifdef FEMCOMPARISON_TIMER
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        solveTimeVec.push_back(static_cast<unsigned long int>(elapsed.count()));
        //solveTimeVec.push_back(static_cast<double>(timer.seconds()));
    }
#endif


    ////PostProcess
    if(pConfig.debugger) {

        ////Calculo do erro
        std::cout << "Computing Error MIXED " << std::endl;

#ifndef OPTMIZE_RUN_TIME
        an.SetExact(config.exact.operator*().ExactSolution());
#endif
        
        std::cout << "DOF = " << cmesh_Mixed->NEquations() << std::endl;

        StockErrors(an,cmesh_Mixed,pConfig.Erro,pConfig.Log,pConfig);

        int dim = config.gmesh->Dimension();
        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        int resolution = 0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.resize(pConfig.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);

    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*pConfig.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(pConfig.h) - log10(pConfig.hLog));
            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << pConfig.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < pConfig.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}

void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.resize(pConfig.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);
    std::cout<<"nnnnnnnn"<<std::endl;
    for(int i=0;i<Errors.size();i++)
        std::cout<<Errors[i]<<std::endl;
    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*pConfig.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(pConfig.h) - log10(pConfig.hLog));
            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << pConfig.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < pConfig.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}
