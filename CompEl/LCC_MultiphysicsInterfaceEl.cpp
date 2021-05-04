//
// Created by victor on 03/05/2021.
//

#include "LCC_MultiphysicsInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "TPZMaterial.h"
#include "pzmultiphysicselement.h"
#include "tpzintpoints.h"

#include "pzmultiphysicscompel.h"
#include "pzgeoel.h"

#include "pzgraphel.h"
#include "pzgraphmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "LCC_LagrangeMultiplier.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("StiffnessInterface");
#endif

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(): TPZMultiphysicsInterfaceElement(), fComputeStiffnessContribution(true)
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right) :
TPZMultiphysicsInterfaceElement(mesh, ref, index,  left, right), fComputeStiffnessContribution(true)
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index):
TPZMultiphysicsInterfaceElement(mesh, ref, index), fComputeStiffnessContribution(true)
{

}

/** @brief create a copy of the given element */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy):
TPZMultiphysicsInterfaceElement(mesh, copy)
{
    fComputeStiffnessContribution = copy.fComputeStiffnessContribution;
}

/** @brief create a copy of the given element using index mapping */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t> & gl2lcElMap):
TPZMultiphysicsInterfaceElement(mesh, copy, gl2lcConMap,gl2lcElMap)
{
    fComputeStiffnessContribution = copy.fComputeStiffnessContribution;
}

void LCC_TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    TPZMaterial *material = this->Material();

    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ek.Reset();
        ef.Reset();
        return;
    }


    InitializeElementMatrix(ek,ef);

    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();
#ifdef PZDEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    std::map<int,TPZMaterialData> datavecleft;
    std::map<int,TPZMaterialData> datavecright;
    TPZMaterialData data;
    InitMaterialData(data, datavecleft, datavecright);

    TPZManVector<TPZTransform<REAL>,6> leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);

    for(auto &it : datavecleft){
        it.second.p = 0;
        int id = it.first;
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(leftel->Element(id));
        if (msp)
        {
            datavecleft[id].p =msp->MaxOrder();
        }
    }
//    for(auto &it : datavecright){
//        it.second.p = 0;
//        int id = it.first;
//        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(rightel->Element(id));
//        if (msp)
//        {
//            datavecright[id].p = msp->MaxOrder();
//        }
//    }

    TPZManVector<int> intleftorder;
    leftel->PolynomialOrder(intleftorder);
    TPZManVector<int> intrightorder;
    rightel->PolynomialOrder(intrightorder);
    int integrationorder = material->GetIntegrationOrder(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);

    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform<> trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    // compute the transformation between neighbours
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);

    TPZTransform<> leftloctr = leftgel->SideToSideTransform(neighleft.Side(), leftgel->NSides()-1);
    TPZTransform<> rightloctr = rightgel->SideToSideTransform(neighright.Side(), rightgel->NSides()-1);
    // transform from the element to the interior of the neighbours
    trleft = leftloctr.Multiply(trleft);
    trright = rightloctr.Multiply(trright);


    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        data.intLocPtIndex = ip;
        intrule->Point(ip, Point, weight);
        ComputeRequiredData(data, Point);
        weight *= fabs(data.detjac);
        trleft.Apply(Point, leftPoint);
        leftel->ComputeRequiredData(leftPoint, leftcomptr, datavecleft);
        trright.Apply(Point, rightPoint);
        rightel->ComputeRequiredData(rightPoint, rightcomptr, datavecright);

        data.x = datavecleft.begin()->second.x;
        material->ContributeInterface(data, datavecleft, datavecright, weight, ek.fMat, ef.fMat);
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nleft element corner nodes: (";
        for(int i = 0 ; i < leftgel->NCornerNodes(); i++){
            sout << leftgel->NodeIndex(i) << " ";
        }
        sout << ")\nright element corner nodes: (";
        for(int i = 0 ; i < rightgel->NCornerNodes(); i++){
            sout << rightgel->NodeIndex(i) << " ";
        }
        sout << ")\n";
        LCC_LagrangeMultiplier *lagMat = dynamic_cast<LCC_LagrangeMultiplier *>(material);
        if(lagMat){
            sout << "fMultiplier = " << lagMat->Multiplier() << "\n";
        }
        ek.fMat.Print(sout);
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

}//CalcStiff
