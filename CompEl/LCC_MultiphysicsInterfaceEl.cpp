//
// Created by victor on 03/05/2021.
//

#include "LCC_MultiphysicsInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(): TPZMultiphysicsInterfaceElement()
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right) :
TPZMultiphysicsInterfaceElement(mesh, ref, index,  left, right)
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index):
TPZMultiphysicsInterfaceElement(mesh, ref, index)
{

}

/** @brief create a copy of the given element */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy):
TPZMultiphysicsInterfaceElement(mesh, copy)
{

}

/** @brief create a copy of the given element using index mapping */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t> & gl2lcElMap):
TPZMultiphysicsInterfaceElement(mesh, copy, gl2lcConMap,gl2lcElMap)
{

}