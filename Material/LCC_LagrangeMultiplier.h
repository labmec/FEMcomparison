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
#include "TPZMaterial.h"
#include "pzvec.h"

/// Material which implements a Lagrange Multiplier
class LCC_LagrangeMultiplier : public TPZMaterial
{

    /** @brief Pointer to a blocked matrix object*/
    TPZManVector<TPZFNMatrix<1000, STATE>,8> fStiffnessMatrix;

    /// Number of state variables
    int fNStateVariables;
    
    /// Dimensiona associated with the material
    int fDimension;
    
    STATE fMultiplier;
    
    public :
	/** @brief Simple constructor */
    LCC_LagrangeMultiplier() : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZMaterial()
    {
        fStiffnessMatrix.Resize(14);
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix].Resize(0,0);
        }
    }
	/** @brief Constructor with the index of the material object within the vector */
    LCC_LagrangeMultiplier(int nummat, int dimension, int nstate) : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZMaterial(nummat), fNStateVariables(nstate), fDimension(dimension), fMultiplier(1.)
    {
        fStiffnessMatrix.Resize(14);
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix].Resize(0,0);
        }
    }
	
	/** @brief Copy constructor */
	LCC_LagrangeMultiplier(const LCC_LagrangeMultiplier &copy) : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZMaterial(copy), fNStateVariables(copy.fNStateVariables), fDimension(copy.fDimension), fMultiplier(copy.fMultiplier)
    {
	    if(fStiffnessMatrix.size() != 14){
	        DebugStop();
	    }
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix] = copy.fStiffnessMatrix[iMatrix];
        }
    }
    
    LCC_LagrangeMultiplier &operator=(const LCC_LagrangeMultiplier &copy)
    {
        TPZMaterial::operator=(copy);
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        if(fStiffnessMatrix.size() != 14){
            DebugStop();
        }
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix] = copy.fStiffnessMatrix[iMatrix];
        }
        return *this;
    }
    
    TPZMaterial *NewMaterial() override
    {
        return new LCC_LagrangeMultiplier(*this);
    }
    
	/** @brief Destructor */
	virtual ~LCC_LagrangeMultiplier()
    {
        
    }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override
    {
        return fDimension;
    }
    
    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }
    
    STATE Multiplier()
    {
        return fMultiplier;
    }
    
    int NStateVariables()
    {
        return fNStateVariables;
    }
	
	virtual std::string Name() override
    {
        return "LCC_LagrangeMultiplier";
    }

    void GetStiffnessMatrix(TPZFMatrix<STATE> &S, int iMatrix){
        if(iMatrix <0 || iMatrix >13){
            DebugStop();
        }
        S = fStiffnessMatrix[iMatrix];
    }

    void FillStiffnessMatrix(const TPZFMatrix<STATE> &S,int iMatrix){
        if(iMatrix <0 || iMatrix >13){
            DebugStop();
        }
        fStiffnessMatrix[iMatrix] = S;
    }
	
    // print the data in human readable form
    virtual void Print(std::ostream &out) override;
	/**
	 * @brief Fill material data parameter with necessary requirements for the ContributeInterface method.
     * @since April 10, 2007
	 */
	/**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavec_left, std::map<int, TPZMaterialData> &datavec_right) override
    {
        data.SetAllRequirements(false);
//        data.fNeedsNormal = true;
        for(auto &it : datavec_left)
        {
            it.second.SetAllRequirements(false);
        }
        for(auto &it : datavec_right)
        {
            it.second.SetAllRequirements(false);
        }        
    }
	
    /**
     * @{
     * @name Contribute methods
     * @}
     */
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override
    {
        DebugStop();
    }
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override
    {
        DebugStop();
    }
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)  override {
        DebugStop();
    }
    
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef)  override {
		DebugStop();
	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		DebugStop();
	}
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    
//    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
//                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
//    {
//        DebugStop();
//    }
	
	/**
	 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
    virtual void ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
	
	/**
	 * @brief It computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;
	
    
	/**
	 * @brief Computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
	virtual void ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef) override
    {
        ContributeInterface(data, dataleft[0], dataright[0], weight, ef);
    }
	
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since February 21, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
    
	
	/**
	 * @brief It computes a contribution to residual vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
	
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
    {
        std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
        DebugStop();
    }
	
    /** @} */
	
	/**
	 * @brief Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
	 * @since Feb 05, 2004
	 */
	/**
	 * @brief Computes interface jump = leftu - rightu
	 * @since Feb 14, 2006
	 */
	virtual void InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump) override
    {
        DebugStop();
    }
	
	
	
	virtual int NStateVariables() const override
    {
        return fNStateVariables;
    }
	
		
	
    /** @{
     * @name Save and Load methods
     */
    
	/** @brief Unique identifier for serialization purposes */
	public:
virtual int ClassId() const override;

	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context) override;
	
    /**
     * @}
     */
};

#endif /* defined(__PZ__LCC_LagrangeMultiplier__) */
