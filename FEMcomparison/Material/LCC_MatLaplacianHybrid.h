#include "TPZMatLaplacianHybrid.h"
class LCC_MatLaplacianHybrid:public TPZMatLaplacianHybrid
{
	public:
    LCC_MatLaplacianHybrid(int matid, int dim);
    
    /*LCC_MatLaplacianHybrid(int matid)
    : TPZRegisterClassId(&LCC_MatLaplacianHybrid::ClassId), LCC_MatLaplacianHybrid(matid)
    {
        
    }*/
    
    LCC_MatLaplacianHybrid();
    
    LCC_MatLaplacianHybrid(const TPZMatLaplacian &copy);
    
    virtual ~LCC_MatLaplacianHybrid();
};
