//
//  LCC_LagrangeMultiplier.h
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#ifndef __PZ__LCC_LagrangeMultiplier__
#define __PZ__LCC_LagrangeMultiplier__

#include <iostream>
#include "TPZLagrangeMultiplierCS.h"
#include "TPZMaterial.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/// Material which implements a Lagrange Multiplier

template<class TVar = STATE>
class LCC_LagrangeMultiplier : public TPZLagrangeMultiplierCS<TVar>
{
public:
    LCC_LagrangeMultiplier() : TPZLagrangeMultiplierCS<TVar>(){

    };
    LCC_LagrangeMultiplier(int nummat, int dimension, int nstate) : TPZLagrangeMultiplierCS<TVar>(nummat,dimension,nstate){

    };
    LCC_LagrangeMultiplier(const LCC_LagrangeMultiplier &copy) : TPZLagrangeMultiplierCS<TVar>(copy){

    };
    /** @brief Destructor */
    virtual ~LCC_LagrangeMultiplier()
    {

    }

    std::string Name() const override
    {return "LCC_LagrangeMultiplierCS";}

    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialDataT<TVar> &data, std::map<int, TPZMaterialData> &datavec_left, std::map<int, TPZMaterialData> &datavec_right) override
    {
        data.SetAllRequirements(false);

        for(auto &it : datavec_left)
        {
            it.second.SetAllRequirements(false);
        }
        for(auto &it : datavec_right)
        {
            it.second.SetAllRequirements(false);
        }
    }

    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                             REAL weight, TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    virtual int ClassId() const override;
};


#endif /* defined(__PZ__LCC_LagrangeMultiplier__) */
