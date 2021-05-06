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
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy):
TPZMultiphysicsInterfaceElement(mesh, copy)
{

}

/** @brief create a copy of the given element using index mapping */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t> & gl2lcElMap):
TPZMultiphysicsInterfaceElement(mesh, copy, gl2lcConMap,gl2lcElMap)
{

}

void LCC_TPZMultiphysicsInterfaceElement::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) {
    TPZMaterial *material = this->Material();

    if (!material) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ek.Reset();
        ef.Reset();
        return;
    }

    LCC_LagrangeMultiplier *lagMat = dynamic_cast<LCC_LagrangeMultiplier *>(material);
    if (!lagMat) {
        DebugStop();
    }

    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();

    InitializeElementMatrix(ek, ef);

    if(this->Dimension() == 1){

        TPZManVector<int> leftNodeIndices;
        int leftNNodes = leftgel->NNodes();
        leftNodeIndices.Resize(leftNNodes,0);

        for (int iNode = 0; iNode < leftgel->NCornerNodes(); iNode++) {
            leftNodeIndices[iNode] = leftgel->NodeIndex(iNode);
        }

        if(leftNodeIndices[0] < leftNodeIndices[1]){ // Matrix A/MinusA
            TPZFMatrix<STATE> A;
            lagMat->GetA(A);
            if(A.Cols() == 0){
                A.Resize(ek.fMat.Rows(),ek.fMat.Cols());
                ComputingCalcStiff(ek,ef);
                lagMat->FillA(ek.fMat);
            }
            else{
#ifdef FEMCOMPARISON_DEBUG
                ComputingCalcStiff(ek,ef);
                if(ek.fMat.Cols() != A.Cols() || ek.fMat.Rows() != A.Rows()){
                    DebugStop();
                }
                STATE relDiff;
                for(int iRow = 0; iRow < A.Rows(); iRow++){
                    for(int iCol = 0 ; iCol < A.Cols(); iCol++){
                        relDiff = (A(iRow,iCol) - ek.fMat(iRow,iCol))/A(iRow,iCol);
                        if(abs(relDiff) > 0.000001){
                            DebugStop();
                        }
                    }
                }
#endif
                ek.fMat = A;
            }
        }
        else{                                        // Matrix B/MinusB
            TPZFMatrix<STATE> B;
            lagMat->GetB(B);
            if(B.Cols() == 0){
                B.Resize(ek.fMat.Rows(),ek.fMat.Cols());
                ComputingCalcStiff(ek,ef);
                lagMat->FillB(ek.fMat);
            }
            else{
#ifdef FEMCOMPARISON_DEBUG
                ComputingCalcStiff(ek,ef);
                if(ek.fMat.Cols() != B.Cols() || ek.fMat.Rows() != B.Rows()){
                    DebugStop();
                }
                STATE relDiff;
                for(int iRow = 0; iRow < B.Rows(); iRow++){
                    for(int iCol = 0 ; iCol < B.Cols(); iCol++){
                        relDiff = (B(iRow,iCol) - ek.fMat(iRow,iCol))/B(iRow,iCol);
                        if(abs(relDiff) > 0.000001){
                            DebugStop();
                        }
                    }
                }
#endif
                ek.fMat = B;
            }
        }
    }
    else{
        ComputingCalcStiff(ek, ef);
    }
}

void LCC_TPZMultiphysicsInterfaceElement::ComputingCalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) {

    TPZMaterial *material = this->Material();

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
        if(leftgel->Dimension() == 1) {
            std::stringstream sout;
            sout << "\nH1-bound element:\n    Element id:" << leftgel->Id() << "\n    Geometric nodes: (";
            for (int i = 0; i < leftgel->NCornerNodes(); i++) {
                sout << leftgel->NodeIndex(i) << " ";
            }
            sout << ")\n    CONNECT INFO. BY GEOMETRIC ELEMENT:\n";
            TPZManVector<REAL, 3> xicenter(gel->Dimension(), 0.);
            TPZManVector<REAL> xcenter(3, 0.);
            for (int i = 0; i < leftgel->NSides(); i++) {
                leftgel->CenterPoint(i, xicenter);
                leftgel->X(xicenter, xcenter);
                sout << "        Cord = [" << xcenter;
                int iCon = i;
                if (leftel->NConnects()) {
                    TPZConnect &con = leftel->Connect(iCon);
                    sout << "] Side = " << i;
                    if (1) {
                        sout << " Seqnumber = " << con.SequenceNumber();
                    }
                    sout << " Order = " << (int) con.Order() << " NState = " << (int) con.NState()
                         << " NShape " << con.NShape();

                    if (1) {
                        sout << " IsLagrangeMult = " << (int) con.LagrangeMultiplier();
                    }
                    sout << " IsCondensed: " << (int) con.IsCondensed();

                    sout << '\n';
                }
            }


            TPZManVector<REAL,3> Node1(leftgel->Dimension()), Node2(leftgel->Dimension());
            Node1[0] = -1.;
            Node2[0] = +1.;

            TPZManVector<REAL,3> X1(3,0), X2(3,0);
            leftgel->X(Node1,X1);
            leftgel->X(Node2,X2);

            leftel->AffineTransform(leftcomptr);
            leftel->ComputeRequiredData(Node1, leftcomptr, datavecleft);

            sout << "    SHAPE FUNCTION VALUE PER NODE:\n";
            sout << "        Node coord:  (" << X1[0] << ", " << X1[1] << ", " << X1[2] << "): \n";
            leftel->AffineTransform(leftcomptr);
            leftel->ComputeRequiredData(Node1, leftcomptr, datavecleft);

            for(auto &it : datavecleft){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }
            sout << "        Node coord:  (" << X2[0] << ", " << X2[1] << ", " << X2[2] << "): \n";
            leftel->AffineTransform(leftcomptr);
            leftel->ComputeRequiredData(Node2, leftcomptr, datavecleft);
            for(auto &it : datavecleft){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }

            sout << "\nH(div)-bound element:\n    Element id:" << rightgel->Id() << "\n    Geometric nodes: (";
            for (int i = 0; i < rightgel->NCornerNodes(); i++) {
                sout << rightgel->NodeIndex(i) << " ";
            }
            sout << ")\n    CONNECT INFO. BY GEOMETRIC ELEMENT:\n";
            int nCon = rightel->NConnects();
            int nSides = gel->NSides();
            int firstSideToHaveConnect = 0;
            if (nSides != nCon) {
                firstSideToHaveConnect = gel->NCornerNodes();
            }
            for (int i = firstSideToHaveConnect; i < nSides; i++) {
                rightgel->CenterPoint(i, xicenter);
                rightgel->X(xicenter, xcenter);
                sout << "        Cord = [" << xcenter;

                int iCon = i - firstSideToHaveConnect;
                if (nCon) {
                    TPZConnect &con = rightel->Connect(iCon);
                    sout << "] Side = " << i;
                    if (1) {
                        sout << " Seqnumber = " << con.SequenceNumber();
                    }
                    sout << " Order = " << (int) con.Order() << " NState = " << (int) con.NState()
                         << " NShape " << con.NShape();

                    if (1) {
                        sout << " IsLagrangeMult = " << (int) con.LagrangeMultiplier();
                    }
                    sout << " IsCondensed: " << (int) con.IsCondensed();

                    sout << '\n';
                }
            }

            sout << "    SHAPE FUNCTION VALUE PER NODE:\n";
            sout << "        Node coord:  (" << X1[0] << ", " << X1[1] << ", " << X1[2] << "): \n";
            rightel->AffineTransform(rightcomptr);
            rightel->ComputeRequiredData(Node1, rightcomptr, datavecright);

            for(auto &it : datavecright){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }
            sout << "        Node coord:  (" << X2[0] << ", " << X2[1] << ", " << X2[2] << "): \n";
            rightel->ComputeRequiredData(Node2, rightcomptr, datavecright);
            for(auto &it : datavecright){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }

            LCC_LagrangeMultiplier *lagMat = dynamic_cast<LCC_LagrangeMultiplier *>(material);
            if (lagMat) {
                sout << "\nfMultiplier = " << lagMat->Multiplier() << "\n";
            }
            ek.fMat.Print(sout);
            sout
                    << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
    }
#endif

}//CalcStiff
