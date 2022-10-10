#ifndef _SPMV570_
#define _SPMV570_
int global_Var;
#endif
#include <iostream>

using namespace std;

namespace SpMV
{
class SparseMatrix
{
private:
    std::size_t _nrows = 0;
    std::size_t _ncols = 0; // _name implies its private variable
    std::size_t _nnz = 0;
protected:
    // data
public:
    SparseMatrix(std::size_t nrows, std::size_t ncols);
    ~SparseMatrix();
};

SparseMatrix::SparseMatrix(const size_t nrows,const std::size_t ncols) : 
    _nrows(nrows), _ncols(ncols)
{
    this->_nnz = 2;
    std::cout << "hiiii!" << std::endl;
    std::cout << "\tnrows: " << this->_nrows << std::endl;
    std::cout << "\tncols: " << this->_ncols << std::endl;
    std::cout << "\tnnz: " << this->_nnz << std::endl;
}
SparseMatrix::~SparseMatrix(/*args*/){
    std::cout << "Destoryah!" << std::endl;
}
};