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
#include <sstream>
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"


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
        TPZNullMaterialCS<STATE> *mat = new TPZNullMaterialCS<STATE>(fWrapId, meshDim - 1,1);
        mcmesh->InsertMaterialObject(mat);
    }
    return;
}

void TPZCreateHybridizedMixedSpaces::AddInterfaceMaterial(TPZMultiphysicsCompMesh *mcmesh){

    int meshDim = mcmesh->Dimension();

    {
        TPZLagrangeMultiplierCS<STATE> *mat = new TPZLagrangeMultiplierCS<STATE>(fInterfaceId, meshDim - 1);
        mcmesh->InsertMaterialObject(mat);
    }
}


void TPZCreateHybridizedMixedSpaces::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mcmesh){
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
        TPZGeoElSide neighLag = gelside.HasNeighbour(fLagrangeId);
        if(!neighLag){
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide cneigh = neighLag.Reference();

        TPZGeoElSide ginterface = gelside.Neighbour();
        if(ginterface.Element()->MaterialId() != fInterfaceId){
            DebugStop();
        }
        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mcmesh,ginterface.Element(),celside,cneigh);
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

    fLagrangeId = matid_base+2*base;
    fWrapId = matid_base;
    fInterfaceId = matid_base + base;
}

void TPZCreateHybridizedMixedSpaces::Print(std::ostream &ofs){
    std::stringstream ss;
    ss << "\n";
    int numStars = 70;
    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";
    ss << __PRETTY_FUNCTION__ << "\nHybridized Mixed approximation space creation set up:\n\n";
    ss << "fSpaceType:\t" <<SpaceTypeName() << "\nNormal Flux Order: " << fNormalFluxOrder;
    ss << "\nPressure Order:\t" << fPressureOrder << "\nfHybridizeBC: " << fHybridizeBC;
    ss << "\n\nmatids:\n\t{";
    for(auto it : fMaterialIds){
        ss << it <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss << "bcmatids:\n\t{";
    for(auto it : fBCMaterialIds){
        ss << it <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss << "peripheralMatids:\n";
    ss << "\twrapid = " << fWrapId;
    ss << "\n\tinterfaceid = " << fInterfaceId;
    ss << "\n\tlagrangeid = " << fLagrangeId <<"\n";

    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    ofs << ss.str();
}

std::string TPZCreateHybridizedMixedSpaces::SpaceTypeName(){
    std::string spaceName;
    switch (fSpaceType) {
        case 0:
            spaceName = "ENone";
            break;
        case 1:
            spaceName = "EHybridizedMixed";
            break;
    }
    return spaceName;
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
        Print(std::cout);
        mcmesh->BuildMultiphysicsSpace(meshvec);

        AddInterfaceMaterial(mcmesh);
        AddInterfaceComputationalElements(mcmesh);

        mcmesh->InitializeBlock();

        std::ofstream ofsmulti("beforeCondensation.txt");
        mcmesh->Print(ofsmulti);
        GroupandCondenseElements(mcmesh);
        std::ofstream ofsmulti2("afterCondensation.txt");
        mcmesh->Print(ofsmulti2);

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

void TPZCreateHybridizedMixedSpaces::GroupandCondenseElements(TPZMultiphysicsCompMesh *cmesh)
{
    /// same procedure as hybridize hdiv
    int64_t nel = cmesh->NElements();
    cmesh->LoadReferences();

    TPZCompEl *cel;
    TPZGeoEl *gel;
    TPZGeoElSide *gelside;
    TPZGeoElSide neigh;
    TPZCompEl *celTarget;

    std::set<int64_t> targetMatIds;
    targetMatIds.insert(fWrapId);
    targetMatIds.insert(fInterfaceId);
    for(auto it : fBCMaterialIds){
        targetMatIds.insert(it);
    }

    std::ofstream ofs1("checkGelmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,ofs1);

    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        cel = cmesh->Element(el);
        if(!cel){
            continue;
        }
        if (cmesh->Element(el)->Dimension() != fDimension) {
            continue;
        }
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);

        gel = cel->Reference();

        int firstSide;
        if (gel->Dimension() >= 0) {
            firstSide = gel->FirstSide(gel->Dimension() - 1);
        } else {
            DebugStop();
        }

        elgr->AddElement(cel);
        for (int iside = firstSide; iside < gel->NSides() - 1; iside++) {
            gelside = new TPZGeoElSide(gel, iside);
            neigh = *gelside;
            if(gelside->HasNeighbour(fBCMaterialIds)){
                neigh = gelside->Neighbour();
                if(fBCMaterialIds.find(neigh.Element()->MaterialId()) != fBCMaterialIds.end()){
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                }
            } else{
                neigh = gelside->Neighbour();
                std::vector<int> targetMatIds= {fWrapId,fInterfaceId};
                for (int ineigh = 0; ineigh < 2; ineigh++) {
                    if (neigh.Element()->MaterialId() != targetMatIds[ineigh]) {
                        DebugStop();
                    }
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                    neigh = neigh.Neighbour();
                }
            }
            gelside = NULL;
        }
    }
    cmesh->ComputeNodElCon();
    if(fSpaceType == EHybridizedMixed)
    {
        int64_t nconnects = cmesh->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == fConfigHybridizedMixed.lagavg) c.IncrementElConnected();
        }
    }
    int neoNel = cmesh->NElements();
    for (int64_t el = nel; el < neoNel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (!elgr) {
            DebugStop();
        }
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
        cond->SetKeepMatrix(false);
        }
}
