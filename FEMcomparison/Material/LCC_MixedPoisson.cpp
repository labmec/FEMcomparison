/**
 * @file
 * @brief Contains the methods of the TPZMixedPoisson class (multiphysics environment)
 * @author Agnaldo Farias
 * @date 2012/05/28
 */

#include "LCC_MixedPoisson.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include <TPZTimer.h>
#ifdef USING_MKL
#include "mkl.h"
#endif

#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif
extern double contributeTime;


LCCMixedPoisson::LCCMixedPoisson(): TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMixedPoisson() {
}

LCCMixedPoisson::LCCMixedPoisson(int matid, int dim): TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMixedPoisson(matid,dim) {
    if (dim < 1) {
        DebugStop();
    }
}

LCCMixedPoisson::~LCCMixedPoisson() {
}

LCCMixedPoisson::LCCMixedPoisson(const TPZMixedPoisson &cp) :TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMixedPoisson(cp) {}

LCCMixedPoisson & LCCMixedPoisson::operator=(const TPZMixedPoisson &copy){
    TPZMixedPoisson::operator=(copy);
    return *this;
}




void LCCMixedPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

    if(fIsStabilized)
        DebugStop();
    TPZTimer timer;
    timer.start();
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    
    STATE force = ff;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    
    TPZFNMatrix<9,STATE> PermTensor;
    TPZFNMatrix<9,STATE> InvPermTensor;
    
    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);
        
    
    if(fUseHdois==true){
        DebugStop();
    }
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
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
#ifdef USING_MKL
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
        alpha = weight;
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
        alpha = weight;
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
            ek(iq,jq) += InvPermTensor(0,0)*weight*prod1;
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
    timer.stop();
    contributeTime += timer.seconds();
}



void LCCMixedPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
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

    REAL v2 = bc.Val2()(0,0);
    REAL v1 = bc.Val1()(0,0);
    REAL u_D = 0;
    REAL normflux = 0.;
    
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
        TPZFNMatrix<9,STATE> PermTensor, InvPermTensor;
        GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);
        
        
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<dim; j++)
            {
                normflux += datavec[0].normal[i]*PermTensor(i,j)*gradu(j,0);
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
        v2 = bc.Val2()(0,0);
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
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
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
                    ef(iq,0)+= gBigNumber*normflux*phiQ(iq,0)*weight;
                    for (int jq=0; jq<phrq; jq++) {
                        
                        ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
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
    
}

int LCCMixedPoisson::ClassId() const{
    return Hash("LCCMixedPoisson") ^ TPZMatPoisson3d::ClassId() << 1;
}

