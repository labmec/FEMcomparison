/**
 * @file
 * @brief Contains the methods of the TPZMixedDarcyFlow class (multiphysics environment)
 * @author Agnaldo Farias
 * @date 2012/05/28
 */

#include "LCC_MixedPoisson.h"
#include "pzlog.h"
#include "TPZBndCondT.h"
#include "TPZMaterialDataT.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#ifdef FEMCOMPARISON_USING_MKL
#include "mkl.h"
#endif
#include "TPZTimer.h"
#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif

#ifdef PZ_LOG
static TPZLogger loggerCTM("contributeTimeVol");
static TPZLogger loggerCTB("contributeTimeBoundary");
#endif

LCCMixedPoisson::LCCMixedPoisson(): TPZRegisterClassId(&LCCMixedPoisson::ClassId), TPZMixedDarcyFlow() {
}

LCCMixedPoisson::LCCMixedPoisson(int matid, int dim): TPZRegisterClassId(&LCCMixedPoisson::ClassId), TPZMixedDarcyFlow(matid,dim) {
    if (dim < 1) {
        DebugStop();
    }
}

LCCMixedPoisson::~LCCMixedPoisson() {
}

LCCMixedPoisson::LCCMixedPoisson(const LCCMixedPoisson &cp) :TPZRegisterClassId(&LCCMixedPoisson::ClassId), TPZMixedDarcyFlow(cp) {}

LCCMixedPoisson & LCCMixedPoisson::operator=(const LCCMixedPoisson &copy){
    TPZMixedDarcyFlow::operator=(copy);
    return *this;
}


void LCCMixedPoisson::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
#ifdef FEMCOMPARISON_TIMER
    extern double contributeTimeVol;
    extern int64_t contributeMaterialCounter;
    //double start = clock();
#endif
#ifdef PZ_LOG
    TPZTimer timer;
    if(loggerCTM.isInfoEnabled())
    timer.start();
#endif

    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    
    STATE force = 0;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x,res);
        force = res[0];
    }
    
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;

    TPZFNMatrix<9,STATE> PermTensor(3,3);
    TPZFNMatrix<9,STATE> InvPermTensor(3,3);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(inv_perm);

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef FEMCOMPARISON_DEBUG
    if(nactive == 4)
    {
        int phrgb = datavec[2].phi.Rows();
        int phrub = datavec[3].phi.Rows();
        if(phrp+phrq+phrgb+phrub != ek.Rows())
        {
            DebugStop();
        }
    }else
    {
        if(phrp+phrq != ek.Rows())
        {
            DebugStop();
        }
    }
#endif
    //Calculate the matrix contribution for flux. Matrix A
    TPZFNMatrix<3,REAL> ivec(3,phrq,0.);
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        for(int id=0; id<3; id++){
            ivec(id,iq) = datavec[0].fDeformedDirections(id,ivecind)*phiQ(ishapeind,0);
        }
    }
#ifdef FEMCOMPARISON_USING_MKL
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = phrq;
        n = phrq;
        k = 3;
        alpha = weight*InvPermTensor(0,0);
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = 3;
        LDB = 3;
        LDC = ek.Rows();
        A = &ivec(0,0);
        B = &ivec(0,0);
        C = &ek(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    //contribucion matrix B
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = phrq;
        n = phrp;
        k = 1;
        alpha = (-1.)*weight;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = phrq;
        LDB = phrp;
        LDC = ek.Rows();
        A = &divQ(0,0);
        B = &phip(0,0);
        C = &ek(0,phrq);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    //contribucion matrix B^t
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = phrp;
        n = phrq;
        k = 1;
        alpha = (-1.)*weight;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = phrp;
        LDB = phrq;
        LDC = ek.Rows();
        A = &phip(0,0);
        B = &divQ(0,0);
        C = &ek(phrq,0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                       m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
#else
    for(int iq=0; iq<phrq; iq++)
    {
        for (int jq=0; jq<phrq; jq++)
        {
            REAL prod1 = ivec(0,iq)*ivec(0,jq) + ivec(1,iq)*ivec(1,jq) + ivec(2,iq)*ivec(2,jq);
            ek(iq,jq) += inv_perm*weight*prod1;
        }
    }
    
    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divQ(iq,0);
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
        }
    }
#endif
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
    if(nactive == 4)
    {
        for(int ip=0; ip<phrp; ip++)
        {
            ek(phrq+ip,phrq+phrp) += phip(ip,0)*weight;
            ek(phrq+phrp,phrq+ip) += phip(ip,0)*weight;
        }
        ek(phrp+phrq+1,phrq+phrp) += -weight;
        ek(phrq+phrp,phrp+phrq+1) += -weight;
    }
#ifdef PZ_LOG
    if(loggerCTM.isDebugEnabled()){
    timer.stop();
#ifdef FEMCOMPARISON_TIMER
        contributeTimeVol += timer.seconds();
        contributeMaterialCounter++;
#endif
    }
#endif
}


void LCCMixedPoisson::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc){
    
#ifdef FEMCOMPARISON_TIMER
    extern double contributeTimeBoundary;
    extern int64_t contributeBoundaryCounter;
#endif
#ifdef PZ_LOG
    TPZTimer timer;
    if(loggerCTB.isDebugEnabled())
    timer.start();
#endif
#ifdef FEMCOMPARISON_DEBUG
    int nref =  datavec.size();
//    if (nref != 2 ) {
//        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
//        DebugStop();
//    }
//    if (bc.Type() > 2 ) {
//        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
//        DebugStop();
//    }
#endif
    
    
    int dim = Dimension();
    
    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0,0);
    REAL u_D = 0;
    REAL normflux = 0.;
    
    if(bc.HasForcingFunctionBC())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunctionBC()(datavec[0].x,res,gradu);
        TPZFNMatrix<9,STATE> PermTensor, InvPermTensor;
        REAL perm;
                perm = GetPermeability(datavec[0].x);
                PermTensor.Diagonal(perm);
                InvPermTensor.Diagonal(1./perm);
        
        for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim; j++)
            {
                normflux += datavec[0].normal[i]*perm*gradu(j,0);
            }
        }
        
        
        if(bc.Type() == 0||bc.Type() == 4)
        {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        }
        else if(bc.Type() == 1 || bc.Type() == 2)
        {
            v2 = -normflux;
            if(bc.Type() ==2)
            {
                v2 = -res[0]+v2/v1;
            }
        }
        else
        {
            DebugStop();
        }
    }else
    {
        v2 = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :            // Neumann condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= fBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= fBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
        
        case 2 :            // mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
                ef(iq,0) += v2*phiQ(iq,0)*weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq,jq) += weight/v1*phiQ(iq,0)*phiQ(jq,0);
                }
            }
            break;
            
        case 4:
            //this case implemented the general Robin boundary condition
            // sigma.n = Km(u-u_D)+g
            //val1(0,0) = Km
            //val2(1,0) = g
            if(IsZero(bc.Val1()(0,0))){
                
                for(int iq=0; iq<phrq; iq++)
                {
                    ef(iq,0)+= fBigNumber*normflux*phiQ(iq,0)*weight;
                    for (int jq=0; jq<phrq; jq++) {
                        
                        ek(iq,jq)+= fBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                    }
                }
                
            }
            else{
            
            REAL InvKm = 1./bc.Val1()(0,0);
            REAL g = normflux;
            for(int in = 0 ; in < phiQ.Rows(); in++) {
                //<(InvKm g - u_D)*(v.n)
                ef(in, 0) +=  (STATE)(InvKm*g -u_D)*phiQ(in,0)* weight;
                for (int jn = 0 ; jn < phiQ.Rows(); jn++) {
                    //InvKm(sigma.n)(v.n)
                    ek(in,jn) += (STATE)(InvKm*phiQ(in,0) * phiQ(jn,0) * weight);
                }
            }
            }
        
            break;
        
    }
#ifdef PZ_LOG
    timer.stop();
#ifdef FEMCOMPARISON_TIMER
    if(loggerCTB.isDebugEnabled()){
    contributeTimeBoundary += timer.seconds();
    contributeBoundaryCounter++;
    }
#endif
#endif
}

int LCCMixedPoisson::ClassId() const{
    return Hash("LCCMixedPoisson") ^ TPZMixedDarcyFlow::ClassId() << 1;
}

