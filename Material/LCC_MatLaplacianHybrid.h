//
//  LCC_MatLaplacianHybrid.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#ifndef TPZMatLaplacianHybrid_hpp
#define TPZMatLaplacianHybrid_hpp

#include <stdio.h>
#include "DarcyFlow/TPZHybridDarcyFlow.h"

class TPZMatLaplacianHybrid : public TPZHybridDarcyFlow
{
    
public:
    
    TPZMatLaplacianHybrid(int matid, int dim);
    
        

    
    TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid(const TPZMatLaplacianHybrid &copy);
    
    virtual ~TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid &operator=(const TPZMatLaplacianHybrid &copy);
    
    virtual TPZMaterial *NewMaterial() override;
    
    virtual int VariableIndex(const std::string &name) override;
    int NSolutionVariables(int var)override;
    
    virtual int NEvalErrors()  override {return 4;}
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)override;

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors)override;

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc) override{

    }

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors,TPZBndCond &bc) override{

    }

    virtual int ClassId() const override;

    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;

    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for (int iref = 0; iref <nref; iref++) {
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsSol = false;
        }
        datavec[0].fNeedsNormal = true;
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = false;
            }
        }
    }
    
    
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    virtual void Read(TPZStream &buf, void *context) override;
    

    
};

#endif /* TPZMatLaplacianHybrid_hpp */
