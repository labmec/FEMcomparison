//
//  LCC_MatLaplacianHybrid.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#ifndef LCC_MatLaplacianHybrid_hpp
#define LCC_MatLaplacianHybrid_hpp

#include <stdio.h>
#include "DarcyFlow/TPZHybridDarcyFlow.h"

class LCC_MatLaplacianHybrid : public TPZHybridDarcyFlow
{
    
public:
    
    LCC_MatLaplacianHybrid(int matid, int dim);
    
    LCC_MatLaplacianHybrid();
    
    LCC_MatLaplacianHybrid(const LCC_MatLaplacianHybrid &copy);
    
    virtual ~LCC_MatLaplacianHybrid();
    
    LCC_MatLaplacianHybrid &operator=(const LCC_MatLaplacianHybrid &copy);
    
    virtual TPZMaterial *NewMaterial() const override;
    
    virtual int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    
    virtual int NEvalErrors() const override {return 4;}
    
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ef) override;

    
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)override;

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc) {

    }

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors,TPZBndCond &bc) {

    }

    virtual int ClassId() const override;

    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsSol = false;
        }
        datavec[0].fNeedsNormal = true;
    }
    
    
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    virtual void Read(TPZStream &buf, void *context) override;
    

    
};

#endif /* LCC_MatLaplacianHybrid_hpp */
