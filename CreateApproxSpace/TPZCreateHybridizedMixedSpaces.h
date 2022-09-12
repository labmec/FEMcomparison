//
// Created by victor on 12/09/2022.
//

#ifndef FEMCOMPARISON_TPZCREATEHYBRIDIZEDMIXEDSPACES_H
#define FEMCOMPARISON_TPZCREATEHYBRIDIZEDMIXEDSPACES_H

#include "TPZMultiphysicsCompMesh.h"
#include <fstream>
#include <ios>

class TPZCreateHybridizedMixedSpaces {
public:
    /// types of spaces this class can create
    enum MSpaceType {ENone, EHybridizedMixed};

private:
    /// the type of space this object will generate
    MSpaceType fSpaceType = EHybridizedMixed;

    int fWrapId;

    int fInterfaceId;

    int fLagrangeId;

    /// the materialids which will be used to create the atomic meshes
    std::set<int> fMaterialIds;

    /// the boundary condition material ids
    std::set<int> fBCMaterialIds;

    /// default internal order for the normal component of HDiv functions
    int fNormalFluxOrder = 1;

    /// default Lagrange multiplier order
    int fPressureOrder = 1;

    /// the dimension of the geometric elements that will be used to generate computational elements
    int fDimension = -1;

    /// the geometric mesh which will generate the computational mesh
    TPZGeoMesh *fGeoMesh = 0;

    /// indicated whether the boundary conditions should be hybridized as well
    bool fHybridizeBC = 0;

    void SetPeripheralsids( int wrapId,int interfaceid,int lagrangeId);

    void DefaultConstructor(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds);

    void ComputePeriferalMaterialIds(int base =10);

public:

    /// All parameters needed for creating a hybrid H1 space
    struct TConfigHybridizedMixed {
        char lagHDiv = 0;
        char lagL2 = 1;
        char lagg = 2;
        char lagavg = 3;

        /// default constructor
        TConfigHybridizedMixed(){}

        /// copy constructor
        TConfigHybridizedMixed(const TConfigHybridizedMixed &copy) = default;

        /// copy operator
        TConfigHybridizedMixed &operator=(const TConfigHybridizedMixed &copy) = default;
    };

    TConfigHybridizedMixed fConfigHybridizedMixed;

    TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds);

    TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &fBCmatIds, int &wrapId,int &interfaceid,int &lagrangeId);

    void SetPOrder(int order){
        fPressureOrder = order;
    }

    void SetNormalFluxOrder(int order){
        fNormalFluxOrder = order;
    }

    TPZMultiphysicsCompMesh* GenerateMesh();

    void Print(std::ostream &out = std::cout);

protected:
    void ConditioningGeomesh();

    TPZCompMesh* CreateDiscHDivMesh();

    TPZCompMesh* CreateL2Mesh();

    TPZCompMesh* CreateConstantMesh(const int &lagNum);

    void AddMaterials(TPZMultiphysicsCompMesh *mcmesh);

    void AddPeripherals(TPZMultiphysicsCompMesh *mcmesh);

    std::string SpaceTypeName();
};






#endif //FEMCOMPARISON_TPZCREATEHYBRIDIZEDMIXEDSPACES_H
