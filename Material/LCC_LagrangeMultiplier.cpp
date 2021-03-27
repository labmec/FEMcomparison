//
//  LCC_LagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "LCC_LagrangeMultiplier.h"
#include "pzaxestools.h"
#ifdef USING_MKL
#include "mkl.h"
#endif
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("LagrangeMultipliersData"));
static LoggerPtr logerror(Logger::getLogger("LagrangeMultipliersError"));
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

//Contribution of skeletal elements.
void LCC_LagrangeMultiplier::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int nmesh = datavec.size();
    if (nmesh!=2) DebugStop();

    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();
    
//------- Block of matrix B ------
    int iq, jp;
	for(iq = 0; iq<phrq; iq++) {
		for(jp=0; jp<phrp; jp++) {
            ek(iq, phrq+jp) += fMultiplier*weight*phiQ(iq,0)*phiP(jp,0);
		}
	}
    
    
//------- Block of matrix B^T ------
    int ip, jq;
	for(ip=0; ip<phrp; ip++) {
		for(jq=0; jq<phrq; jq++) {
			ek(ip + phrq,jq) += fMultiplier*weight*phiP(ip,0)*phiQ(jq,0);
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
void LCC_LagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> *phiLPtr = 0, *phiRPtr = 0;
    for (int i=0; i<dataleft.size(); i++) {
        if (dataleft[i].phi.Rows() != 0) {
            phiLPtr = &dataleft[i].phi;
            break;
        }
    }
    for (int i=0; i<dataright.size(); i++) {
        if (dataright[i].phi.Rows() != 0) {
            phiRPtr = &dataright[i].phi;
            break;
        }
    }
    
    if(!phiLPtr || !phiRPtr)
    {
        DebugStop();
    }
    TPZFMatrix<REAL> &phiL = *phiLPtr;
    TPZFMatrix<REAL> &phiR = *phiRPtr;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    static int count  = 0;

    if((nrowl+nrowr)*fNStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nrowl <<
        " nrowr " << nrowr << " may give wrong result " << std::endl;
        count++;
    }

    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
    int il,jl,ir,jr;

#ifdef USING_MKL
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k,phrL,phrR;
        phrL = phiL.Rows();
        phrR = phiR.Rows();
        m = phrL;
        n = phrR;
        k = 1;
        alpha = weight * fMultiplier ;
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
        alpha = weight * fMultiplier ;
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
void LCC_LagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
//	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
//	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	
//	TPZFNMatrix<660> dphiL, dphiR;
//	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
//	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	

	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
#ifdef PZDEBUG
    if(phiL.Rows()*fNStateVariables+phiR.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
	int il,jl,ir,jr;
    
    
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
    
}

/**
 * @brief It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void LCC_LagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

// print the data in human readable form
void LCC_LagrangeMultiplier::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
}

