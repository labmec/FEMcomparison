//
//  LCC_LagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "LCC_LagrangeMultiplier.h"
#include "pzaxestools.h"
#ifdef FEMCOMPARISON_USING_MKL
#include "mkl.h"
#endif
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("LagrangeMultipliersData"));
static LoggerPtr logerror(Logger::getLogger("LagrangeMultipliersError"));
#endif
#include "TPZTimer.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger loggerCTI("contributeTimeInterface");
#endif

/** @brief Unique identifier for serialization purposes */
template<class TVar>
int LCC_LagrangeMultiplier<TVar>::ClassId() const{
    return Hash("LCC_LagrangeMultiplier") ^ TPZMaterial::ClassId() << 1;
}

//Contribution of skeletal elements.
template<class TVar>
void LCC_LagrangeMultiplier<TVar>::Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef)
{
    int nmesh = datavec.size();
    if (nmesh!=2) DebugStop();

    TVar Multiplier = this->Multiplier();

    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();
    
//------- Block of matrix B ------
    int iq, jp;
	for(iq = 0; iq<phrq; iq++) {
		for(jp=0; jp<phrp; jp++) {
            ek(iq, phrq+jp) += Multiplier*weight*phiQ(iq,0)*phiP(jp,0);
		}
	}
    
    
//------- Block of matrix B^T ------
    int ip, jq;
	for(ip=0; ip<phrp; ip++) {
		for(jq=0; jq<phrq; jq++) {
			ek(ip + phrq,jq) += Multiplier*weight*phiP(ip,0)*phiQ(jq,0);
		}
	}
}

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
#ifdef FEMCOMPARISON_TIMER
extern double contributeTimeInterface;
#endif

template<class TVar>
void LCC_LagrangeMultiplier<TVar>::ContributeInterface(const TPZMaterialDataT<TVar> &data,
                         const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                         const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                         REAL weight, TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef)
{
#ifdef PZ_LOG
    TPZTimer timer;
    if(loggerCTI.isDebugEnabled()){
    timer.start();}
#endif
#ifdef FEMCOMPARISON_DEBUG
    if(dataleft.size() != 1 || dataright.size() != 1) DebugStop();
#endif
    TPZFMatrix<REAL> &phiL = dataleft.begin()->second.phi;
    TPZFMatrix<REAL> &phiR = dataright.begin()->second.phi;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    static int count  = 0;

    int NStateVariables = this->NStateVariables();
    TVar Multiplier = this->Multiplier();

    if((nrowl+nrowr)*NStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nrowl <<
        " nrowr " << nrowr << " may give wrong result " << std::endl;
        count++;
    }

    int secondblock = ek.Rows()-phiR.Rows()*NStateVariables;
    int il,jl,ir,jr;

#ifdef FEMCOMPARISON_USING_MKL
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k,phrL,phrR;
        phrL = phiL.Rows();
        phrR = phiR.Rows();
        m = phrL;
        n = phrR;
        k = 1;
        alpha = weight * Multiplier ;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = phrL+phrR;
        LDA = phiL.Rows();
        LDB = 1;
        C = &ek(0,phrL);
        A = &phiL(0,0);
        B = &phiR(0,0);
        //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    {
      double *A, *B, *C;
        double alpha, beta;
        int m,n,k,phrL,phrR;
        phrL = phiL.Rows();
        phrR = phiR.Rows();
        m = phrR;
        n = phrL;
        k = 1;
        alpha = weight * Multiplier ;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = phrL+phrR;
        LDA = phiR.Rows();
        LDB = 1;
        C = &ek(phrL,0);
        B = &phiL(0,0);
        A = &phiR(0,0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    
#else

    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
        }
    }
    
    //	// 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
        }
    }
#endif
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream valuenn;
        ek.Print("ek = ",valuenn,EMathematicaInput);
        ef.Print("ef = ",valuenn,EMathematicaInput);
        LOGPZ_DEBUG(logdata,valuenn.str());
    }
#endif
#ifdef PZ_LOG
    if(loggerCTI.isDebugEnabled())
    timer.stop();
    contributeTimeInterface += timer.seconds();
#endif
}