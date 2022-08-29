//
//  LCC_LagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//
#include "fstream"
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
int LCC_LagrangeMultiplier::ClassId() const{
    return Hash("LCC_LagrangeMultiplier") ^ TPZMaterial::ClassId() << 1;
}

/** @brief Saves the element data to a stream */
void LCC_LagrangeMultiplier::Write(TPZStream &buf, int withclassid) const
{
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

/** @brief Reads the element data from a stream */
void LCC_LagrangeMultiplier::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf, context);
    buf.Read(&fNStateVariables);
    
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

void LCC_LagrangeMultiplier::ContributeInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &dataleft, const std::map<int, TPZMaterialDataT<STATE>> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
#ifdef FEMCOMPARISON_DEBUG
    if(dataleft.size() != 1 || dataright.size() != 1) DebugStop();
#endif

    //const TPZFMatrix<REAL> *phiL = &(dataleft.begin()->second.phi);
    TPZFMatrix<REAL> phiLdummy = dataleft.begin()->second.phi;
    TPZFMatrix<REAL> *phiL = &phiLdummy;

    TPZFMatrix<REAL> phiRdummy = dataright.begin()->second.phi;
    TPZFMatrix<REAL> *phiR = &phiRdummy;

    int nrowl = phiL->Rows();
    int nrowr = phiR->Rows();
    static int count  = 0;

    if((nrowl+nrowr)*fNStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nrowl <<
        " nrowr " << nrowr << " may give wrong result "<< "fNStateVariables\t"<< fNStateVariables << "count " << count  << std::endl;
        count++;
    }

    int secondblock = ek.Rows()-phiR->Rows()*fNStateVariables;
    int il,jl,ir,jr;

#ifdef FEMCOMPARISON_USING_MKL
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k,phrL,phrR;
        phrL = phiL->Rows();
        phrR = phiR->Rows();
        m = phrL;
        n = phrR;
        k = 1;
        alpha = weight * fMultiplier ;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = phrL+phrR;
        LDA = phiL->Rows();
        LDB = 1;
        C = &ek(0,phrL);
        A = &phiLdummy(0,0);
        B = &phiRdummy(0,0);
        //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    {
      double *A, *B, *C;
        double alpha, beta;
        int m,n,k,phrL,phrR;
        phrL = phiL->Rows();
        phrR = phiR->Rows();
        m = phrR;
        n = phrL;
        k = 1;
        alpha = weight * fMultiplier ;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = phrL+phrR;
        LDA = phiR->Rows();
        LDB = 1;
        C = &ek(phrL,0);
        B = &phiLdummy(0,0);
        A = &phiRdummy(0,0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    
#else

    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                double L = *(&phiLdummy(0,0)+il);
                double R = *(&phiRdummy(0,0)+jr);
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * L * R;
            }
        }
    }
    
    //	// 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                double L = *(&phiLdummy(0,0)+jl);
                double R = *(&phiRdummy(0,0)+ir);
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (R * L);

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

}

// print the data in human readable form
void LCC_LagrangeMultiplier::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZLagrangeMultiplierCS<STATE>::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
}

