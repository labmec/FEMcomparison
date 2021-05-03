//
// Created by victor on 03/05/2021.
//

#ifndef FEMCOMPARISON_LCC_MULTIPHYSICSINTERFACEEL_H
#define FEMCOMPARISON_LCC_MULTIPHYSICSINTERFACEEL_H

#include "pzcompel.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"

class LCC_TPZMultiphysicsInterfaceElement : public TPZMultiphysicsInterfaceElement{
public:
    /** @brief Default constructor */
    LCC_TPZMultiphysicsInterfaceElement();

    LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);

    LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index);

    LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy);

    LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t> & gl2lcElMap);

    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        return new LCC_TPZMultiphysicsInterfaceElement(mesh,*this);
    }

    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const override
    {
        return new LCC_TPZMultiphysicsInterfaceElement(mesh,*this,gl2lcConMap,gl2lcElMap);
    }

/** @brief Default destructor */
virtual ~LCC_TPZMultiphysicsInterfaceElement(){

}

};


#endif //FEMCOMPARISON_LCC_MULTIPHYSICSINTERFACEEL_H
