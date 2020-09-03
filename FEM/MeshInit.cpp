//
// Created by victor on 02/09/2020.
//

#include "MeshInit.h"
#include "DataStructure.h"
#include <TPZMultiphysicsCompMesh.h>
#include "mixedpoisson.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"
#include "TPZGenGrid2D.h"

void InsertMaterialMixed_MultiK(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig &config, PreConfig &pConfig){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

    int defaultOrder = config.k+config.n;
    cmesh_mixed -> SetDefaultOrder(defaultOrder);
    cmesh_mixed -> SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    SetFExact(mat1,mat2,pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMixedPoisson *material_Q1 = new TPZMixedPoisson(matID_Q1,dim); //Using standard PermealityTensor = Identity.
    TPZMixedPoisson *material_Q2 = new TPZMixedPoisson(matID_Q2,dim);
    material_Q1->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material_Q1->SetForcingFunctionExact(config.exact.operator*().Exact());
    material_Q2->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material_Q2->SetForcingFunctionExact(config.exact.operator*().Exact());

    material_Q1->SetPermeability(pConfig.perm_Q1);
    material_Q2->SetPermeability(pConfig.perm_Q2);

    cmesh_mixed->InsertMaterialObject(material_Q1);
    cmesh_mixed->InsertMaterialObject(material_Q2);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    BCond0_Q1->SetForcingFunction(config.exact.operator*().Exact());

    TPZMaterial *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    BCond0_Q2->SetForcingFunction(config.exact.operator*().Exact());

    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q1->CreateBC(material_Q2, -9, neumann, val1, val2);

    cmesh_mixed->InsertMaterialObject(BCond0_Q1);
    cmesh_mixed->InsertMaterialObject(BCond0_Q2);
    cmesh_mixed->InsertMaterialObject(BCond1_Q1);
    cmesh_mixed->InsertMaterialObject(BCond1_Q2);
}

void InsertMaterialHybrid_MultiK(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig){

    TPZGeoMesh* gmesh = cmesh_H1Hybrid->Reference();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    SetFExact(mat1,mat2,pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMatLaplacianHybrid *material_Q1 = new TPZMatLaplacianHybrid(matID_Q1, dim);
    TPZMatLaplacianHybrid *material_Q2 = new TPZMatLaplacianHybrid(matID_Q2, dim);

    material_Q1->SetPermeability(pConfig.perm_Q1);
    material_Q2->SetPermeability(pConfig.perm_Q2);

    cmesh_H1Hybrid->InsertMaterialObject(material_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(material_Q2);

    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        material_Q1->SetForcingFunction(mat1->ForcingFunction());
        material_Q1->SetForcingFunctionExact(mat1->Exact());

        material_Q2->SetForcingFunction(mat2->ForcingFunction());
        material_Q2->SetForcingFunctionExact(mat2->Exact());
    }

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    TPZMaterial *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0_Q1->SetForcingFunction(mat1->Exact());
        BCond0_Q2->SetForcingFunction(mat2->Exact());
    }
    val2.Zero();
    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q2);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q2);
}

void SetFExact(TLaplaceExample1 *mat1, TLaplaceExample1 *mat2,PreConfig &pConfig){
    switch(pConfig.type) {
        case 0:
            mat1->fExact = TLaplaceExample1::ESinSin;
            mat2->fExact = TLaplaceExample1::ESinSin;
            break;
        case 1:
            mat1->fExact = TLaplaceExample1::EArcTan;
            mat2->fExact = TLaplaceExample1::EArcTan;
            break;
        case 2:
            mat1->fExact = TLaplaceExample1::ESteklovNonConst;
            mat2->fExact = TLaplaceExample1::ESteklovNonConst;
            break;
        default:
            DebugStop();
            break;
    }
    mat1->fSignConvention = 1;
    mat2->fSignConvention = 1;
}

void SetMultiPermeMaterials(TPZGeoMesh* gmesh){

    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    int numEl = elvec.NElements();
    TPZVec<int64_t> nodeInd;
    int numberNodes;
    TPZManVector<REAL,2> elCenterCoord = {0.,0.};

    //Getting center coordinates of the element
    for(int ind = 0; ind < numEl; ind++) {
        TPZGeoEl *gel = elvec[ind];
        gel->GetNodeIndices(nodeInd);
        numberNodes = gel->NNodes();
        elCenterCoord[0] = elCenterCoord[1] = 0.;
        for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
            TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
            elCenterCoord[0] += nodo.Coord(0);
            elCenterCoord[1] += nodo.Coord(1);
        }
        elCenterCoord[0] /= numberNodes;
        elCenterCoord[1] /= numberNodes;

        // Changing the permeability
        //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
        if((elCenterCoord[0] >=0. && elCenterCoord[1] >= 0.) || (elCenterCoord[0] < 0. && elCenterCoord[1] <= 0.)) { // first/third quadrants
            if(gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(2);}
            else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-5);} // Assumes dirichlet boundary condition
        }
        else{
            if(gel->Dimension() == gmesh->Dimension())  { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(3);}
            else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-6);}// Assumes dirichlet boundary condition
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

TPZGeoMesh* CreateGeoMesh_OriginCentered(int nel, TPZVec<int>& bcids) {

    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, -1.), x1(3, 1.);
    x1[2] = x0[2] = 0.;
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);

    //TPZGenGrid2D gen(nx, x0, x1);
    gen.SetRefpatternElements(true);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, bcids[0]);
    gen.SetBC(gmesh, 5, bcids[1]);
    gen.SetBC(gmesh, 6, bcids[2]);
    gen.SetBC(gmesh, 7, bcids[3]);

    gmesh->SetDimension(2);

    return gmesh;
}

void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_flux->SetDefaultOrder(config.k);
    cmesh_flux->SetDimModel(config.gmesh->Dimension());

    cmesh_flux->SetAllCreateFunctionsHDiv();

    TPZNullMaterial *material = new TPZNullMaterial(matID);
    material->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    //(PROPRIETARY)
    TPZMaterial *BCond0 = material->CreateBC(material,-1,dirichlet,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond0);

    TPZMaterial *BCond1 = material->CreateBC(material,-2,neumann,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond1);

    //This part is for multi-K problemas
    if (pConfig.type == 2) BuildFluxMesh_MultiK(cmesh_flux,config);

    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config) {

    int dim = config.gmesh->Dimension();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dirichlet = 0;
    int neumann = 1;

    TPZNullMaterial *material_Q1 = new TPZNullMaterial(matID_Q1);
    TPZNullMaterial *material_Q2 = new TPZNullMaterial(matID_Q2);
    material_Q1->SetDimension(dim);
    material_Q2->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material_Q1);
    cmesh_flux->InsertMaterialObject(material_Q2);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);

    //(PROPRIETARY)
    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    TPZMaterial *BCond0_Q2 = material_Q1->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond0_Q1);
    cmesh_flux->InsertMaterialObject(BCond0_Q2);

    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q2->CreateBC(material_Q2, -9, neumann, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond1_Q1);
    cmesh_flux->InsertMaterialObject(BCond1_Q2);
}

void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config, PreConfig &pConfig){
    int matID = 1;
    int dim = config.gmesh->Dimension();

    int potential_order = config.k + config.n;
    cmesh_p->SetDefaultOrder(potential_order);
    cmesh_p->SetDimModel(dim);

    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial *material = new TPZNullMaterial(matID); material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);

    if (pConfig.type == 2) BuildPotentialMesh_MultiK(cmesh_p, config);

    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();

    TPZNullMaterial *material_Q1 = new TPZNullMaterial(matID_Q1); material_Q1->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q1);

    TPZNullMaterial *material_Q2 = new TPZNullMaterial(matID_Q2); material_Q2->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q2);
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    if(pConfig.type != 2) {
        cmesh_mixed->SetDefaultOrder(config.k + config.n);
        cmesh_mixed->SetDimModel(dim);
        cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

        TPZMixedPoisson *material = new TPZMixedPoisson(matID, dim); //Using standard PermealityTensor = Identity.
        material->SetForcingFunction(config.exact.operator*().ForcingFunction());
        material->SetForcingFunctionExact(config.exact.operator*().Exact());
        cmesh_mixed->InsertMaterialObject(material);

        //Boundary Conditions
        TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);

        TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
        BCond0->SetForcingFunction(config.exact.operator*().Exact());

        TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        cmesh_mixed->InsertMaterialObject(BCond0);
        cmesh_mixed->InsertMaterialObject(BCond1);
    }

    else {
        InsertMaterialMixed_MultiK(cmesh_mixed, config, pConfig);
    }
}

void InsertMaterialHybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    if(pConfig.type != 2) {
        TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);

        material->SetPermeability(1.);

        cmesh_H1Hybrid->InsertMaterialObject(material);
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            material->SetForcingFunction(
                    config.exact.operator*().ForcingFunction());
            material->SetForcingFunctionExact(config.exact.operator*().Exact());
        }
        //    TPZMaterial * mat(material);
        //    cmesh->InsertMaterialObject(mat);

        // Inserts boundary conditions
        TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
        TPZMaterial *BCond0 =
                material->CreateBC(material, -1, dirichlet, val1, val2);
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            BCond0->SetForcingFunction(config.exact.operator*().Exact());
        }
        val2.Zero();
        TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        cmesh_H1Hybrid->InsertMaterialObject(BCond0);
        cmesh_H1Hybrid->InsertMaterialObject(BCond1);
    }
    else {
        InsertMaterialHybrid_MultiK(cmesh_H1Hybrid, config, pConfig);
    }
}

