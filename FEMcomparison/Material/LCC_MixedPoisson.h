

#ifndef LCC_mixedpoisson_h
#define LCC_mixedpoisson_h

#include "mixedpoisson.h"

/**
 * @ingroup material
 * @author Agnaldo Farias
 * @since 5/28/2012
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -(k/visc)*grad(p)  ==> Int{Q.q}dx - (k/visc)*Int{p*div(q)}dx + (k/visc)*Int{pD(q.n)}ds = 0  (Eq. 1)  \f$
 *
 * \f$ div(Q) = f  ==> Int{div(Q)*v}dx = Int{f*v}dx (Eq. 2) \f$
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */


class LCCMixedPoisson : public TPZMixedPoisson {
    
protected:
  
    
public:
    LCCMixedPoisson();
    
    LCCMixedPoisson(int matid, int dim);
    
    virtual ~LCCMixedPoisson();
    
    LCCMixedPoisson(const TPZMixedPoisson &cp);
    
    LCCMixedPoisson &operator=(const TPZMixedPoisson &copy);
    
    virtual TPZMaterial * NewMaterial() override{
        return new LCCMixedPoisson(*this);
    }
  
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
   
    public:
virtual int ClassId() const  override;

};

#endif

