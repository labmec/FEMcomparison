//
//  LCC_MatLaplacianHybrid.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#ifndef LCC_MatLaplacianHybrid_hpp
#define LCC_MatLaplacianHybrid_hpp

#include <stdio.h>
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


typedef TPZDarcyFlow TPZMatLaplacian;

class LCC_MatLaplacianHybrid  : public TPZMatCombinedSpacesT<STATE>, public TPZMatErrorCombinedSpaces<STATE>, public TPZDarcyFlow
{

public:

    LCC_MatLaplacianHybrid (int matid, int dim);

    LCC_MatLaplacianHybrid ();

    LCC_MatLaplacianHybrid (const TPZMatLaplacian &copy);

    virtual ~LCC_MatLaplacianHybrid ();

    LCC_MatLaplacianHybrid  &operator=(const LCC_MatLaplacianHybrid  &copy);

    virtual TPZMaterial *NewMaterial() const override;

    virtual int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;

    virtual int NEvalErrors() const override {return 4;}

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE> > &datavec, int var, TPZVec<STATE> &Solout)override;

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors) override;

    void VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary);

    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) override;

    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override
    {
        datavec[0].fNeedsNormal = true;

        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsSol = false;
        }

        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = false;
            }
        }
    }


    virtual void ErrorsBC(TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc){
    }

    virtual void ErrorsBC(TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc)
    {
    }


    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    ///  HDiv simulations use an additional integration order.
    ///  In order to test if HybridH1 and HDiv are equivalents, both integration orders shall be equivalent.
    ///  virtual int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const;

    virtual int ClassId() const override;


    virtual void Write(TPZStream &buf, int withclassid) const override;

    virtual void Read(TPZStream &buf, void *context) override;

    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
                                         int id, int type,
                                         const TPZFMatrix<STATE> &val1,
                                         const TPZVec<STATE> &val2) override
    {
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE> >
                (reference,id, type,val1,val2);
    }


};
#endif /* LCC_MatLaplacianHybrid _hpp */
