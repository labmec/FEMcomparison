//
// Created by victor on 12/09/2022.
//

#include "TPZCreateHybridizedMixedSpaces.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzgeoelside.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"


TPZCreateHybridizedMixedSpaces::TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds){
DefaultConstructor(gmesh,matids,bcmatIds);
ComputePeriferalMaterialIds();
}

TPZCreateHybridizedMixedSpaces::TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &fBCmatIds, int &wrapId,int &interfaceid,int &lagrangeId){
SetPeripheralsids(wrapId,interfaceid,lagrangeId);
DefaultConstructor(gmesh,matids,fBCmatIds);
}

void TPZCreateHybridizedMixedSpaces::SetPeripheralsids( int wrapId,int interfaceId,int lagrangeId){
    fWrapId = wrapId;
    fInterfaceId = interfaceId;
    fLagrangeId = lagrangeId;
}

void TPZCreateHybridizedMixedSpaces::DefaultConstructor(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds){
    fMaterialIds = matids;
    fBCMaterialIds = bcmatIds;
    fGeoMesh = gmesh;
    fDimension = gmesh->Dimension();
}

void TPZCreateHybridizedMixedSpaces::ConditioningGeomesh(){
    int64_t numelem = fGeoMesh->NElements();
    int meshDim = fGeoMesh->Dimension();
    for(int iel =0 ; iel < numelem; iel++){
        TPZGeoEl* gel = fGeoMesh->Element(iel);
        if(!gel || gel->HasSubElement() || gel->Dimension() != meshDim){
            continue;
        }
        int matId = gel->MaterialId();
        if(fMaterialIds.find(matId) == fMaterialIds.end()){
            continue;
        }
        int firstSide = gel->FirstSide(meshDim - 1);
        for(int iside = firstSide; iside < gel->NSides()-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            if(gelside.HasNeighbour(fBCMaterialIds)){
                continue;
            }
            if(!gelside.HasNeighbour(fLagrangeId)){
                TPZGeoElBC gelbc(gelside,fLagrangeId);
            }
            TPZGeoElBC gelbcInt(gelside,fInterfaceId);
            TPZGeoElBC gelbcWrap(gelside,fWrapId);
#ifdef PZDEBUG
            {
                auto neighbour = gelside.Neighbour();
                if(neighbour.Element()->MaterialId() !=fWrapId){
                    DebugStop();
                }
                auto neighneigh = neighbour.Neighbour();
                if(neighneigh.Element()->MaterialId() != fInterfaceId){
                    DebugStop();
                }
                auto neighneighneigh = neighneigh.HasNeighbour(fLagrangeId);
                if(!neighneighneigh){
                    DebugStop();
                }
                if(neighneighneigh.Element()->MaterialId() != fLagrangeId){
                    DebugStop();
                }
            };
#endif
        }
    }
    {// Printing mesh
#ifdef PZDEBUG
        std::ofstream ofs1("gelCreation.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,ofs1);
        return;
#endif
    }

}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateDiscHDivMesh(){
    TPZCompMesh *hdivMesh = new TPZCompMesh(fGeoMesh);
    hdivMesh->SetAllCreateFunctionsHDiv();
    hdivMesh->ApproxSpace().CreateDisconnectedElements(true);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        hdivMesh->InsertMaterialObject(mat);
    }

    for(auto it : fBCMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim - 1);
        hdivMesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(fWrapId, meshDim - 1);
        hdivMesh->InsertMaterialObject(mat);
    }

    hdivMesh->AutoBuild(fMaterialIds);
    int numEl = hdivMesh->NElements();
    for(int iel = 0; iel < numEl; iel++){
        auto *cel = hdivMesh->Element(iel);
        auto *gel = cel->Reference();
        if(!cel){
            DebugStop();
        }
        cel->LoadElementReference();
        int firstSide = gel->FirstSide(meshDim - 1);
        for(int iside = firstSide; iside < gel->NSides()-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            auto neigh = gelside.HasNeighbour(fBCMaterialIds);
            if(neigh){
                hdivMesh->ApproxSpace().CreateCompEl(neigh.Element(),*hdivMesh);
                neigh.Element()->ResetReference();
            }
            else{
                neigh = gelside.Neighbour();
                hdivMesh->ApproxSpace().CreateCompEl(neigh.Element(),*hdivMesh);
                neigh.Element()->ResetReference();
            }
        }
        gel->ResetReference();
    }

    int64_t nc = hdivMesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        hdivMesh->ConnectVec()[icon].SetLagrangeMultiplier(fConfigHybridizedMixed.lagHDiv);
    }
    return hdivMesh;
}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateL2Mesh(){
    TPZCompMesh *L2Mesh = new TPZCompMesh(fGeoMesh);
    L2Mesh->SetAllCreateFunctionsContinuous();
    L2Mesh->ApproxSpace().CreateDisconnectedElements(true);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        L2Mesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(fLagrangeId, meshDim - 1);
        L2Mesh->InsertMaterialObject(mat);
    }

    L2Mesh->AutoBuild(fMaterialIds);

    L2Mesh->SetDimModel(meshDim-1);
    std::set<int> ids = {fLagrangeId};
    L2Mesh->AutoBuild(ids);

    int64_t nc = L2Mesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        L2Mesh->ConnectVec()[icon].SetLagrangeMultiplier(fConfigHybridizedMixed.lagL2);
    }
    return L2Mesh;
}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateConstantMesh(const int &lagNum){
    TPZCompMesh *constMesh = new TPZCompMesh(fGeoMesh);
    constMesh->SetAllCreateFunctionsDiscontinuous();
    constMesh->SetDefaultOrder(0);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        constMesh->InsertMaterialObject(mat);
    }

    constMesh->AutoBuild(fMaterialIds);

    int64_t nc = constMesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        constMesh->ConnectVec()[icon].SetLagrangeMultiplier(lagNum);
    }
    return constMesh;
}

void TPZCreateHybridizedMixedSpaces::AddMaterials(TPZMultiphysicsCompMesh *mcmesh){

    int meshDim = mcmesh->Dimension();

    TPZMixedDarcyFlow* refmat = 0;
    for(auto it : fMaterialIds){
        TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(it,meshDim);
        mcmesh->InsertMaterialObject(mat);
        refmat = mat;
    }

    for(auto it : fBCMaterialIds){
        TPZFMatrix<STATE> val1(1,1,0);
        TPZVec<STATE> val2(1,0);
        auto *mat = refmat->CreateBC(refmat,it,0,val1,val2);
        mcmesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterialCS<STATE> *mat = new TPZNullMaterialCS<STATE>(fLagrangeId, meshDim - 1,1);
        mcmesh->InsertMaterialObject(mat);
    }

    {
        TPZLagrangeMultiplierCS<STATE> *mat = new TPZLagrangeMultiplierCS<STATE>(fInterfaceId, meshDim - 1);
        mcmesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterialCS<STATE> *mat = new TPZNullMaterialCS<STATE>(fWrapId, meshDim - 1,1);
        mcmesh->InsertMaterialObject(mat);
    }
    return;
}

void TPZCreateHybridizedMixedSpaces::AddPeripherals(TPZMultiphysicsCompMesh *mcmesh){
    int numEl = mcmesh->NElements();

    auto fGeoMesh = mcmesh->Reference();
    fGeoMesh->ResetReference();
    mcmesh->LoadReferences();

    for(int iel = 0; iel < numEl; iel++){
        TPZCompEl *cel = mcmesh->Element(iel);
        if(!cel || !cel->Reference()){
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fWrapId){
            continue;
        }
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neigh = gelside.HasNeighbour(fLagrangeId);
        if(!neigh){
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide cneigh = neigh.Reference();

        TPZGeoElSide interside = gelside.HasNeighbour(fInterfaceId);
        if(!interside){
            DebugStop();
        }
        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mcmesh,interside.Element(),celside,cneigh);
    }
}

void TPZCreateHybridizedMixedSpaces::ComputePeriferalMaterialIds(int base)
{
    if(base < 2) base = 2;
    int max_matid = 0;
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) continue;
        max_matid = std::max(max_matid,gel->MaterialId());
    }
    int remain = max_matid % base;
    int matid_base = max_matid-remain + base;

    fLagrangeId = matid_base;
    fWrapId = matid_base+base;
    fInterfaceId = matid_base + 2*base;
}

TPZMultiphysicsCompMesh* TPZCreateHybridizedMixedSpaces::GenerateMesh(){
    {
        ConditioningGeomesh();
        TPZCompMesh *hdivMesh = CreateDiscHDivMesh();
        TPZCompMesh *L2Mesh = CreateL2Mesh();
        TPZCompMesh *gspace   =   CreateConstantMesh(fConfigHybridizedMixed.lagg);
        TPZCompMesh *avgspace = CreateConstantMesh(fConfigHybridizedMixed.lagavg);

        TPZVec<TPZCompMesh*> meshvec = {hdivMesh,L2Mesh,gspace,avgspace};
        TPZMultiphysicsCompMesh *mcmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
        mcmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

        AddMaterials(mcmesh);
        mcmesh->BuildMultiphysicsSpace(meshvec);

        AddPeripherals(mcmesh);

        mcmesh->LoadReferences();
        mcmesh->InitializeBlock();

        {
            std::ofstream ofs1("hdiv.txt"),ofs2("L2.txt"),ofs3("g.txt"),ofs4("avg.txt");
            std::ofstream ofsmulti("multi.txt");

            hdivMesh->Print(ofs1);
            L2Mesh->Print(ofs2);
            gspace->Print(ofs3);
            avgspace->Print(ofs4);

            mcmesh->Print(ofsmulti);
        }

        return mcmesh;
        //Agroupar e condensar
        //percorrer vizinhos, colocar em grupo e condensar
        //Agroupar wrap com interface
        // teremos malha hdiv hibridizada
        // Opção de hibridizar contorno
    }
}
