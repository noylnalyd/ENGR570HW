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

    void setCoefficient(const size_t rowIdx, const size_t colIdx, const double Aij);
    void assembleStorage();
    SparseMatrix getFormat();
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
void SparseMatrix::setCoefficient(const size_t rowIdx, const size_t colIdx, const double Aij){
    cout << "howdy u set a variable <3" << endl;
}
void SparseMatrix::assembleStorage(){
    cout << "assemblin the storage toady?" << endl;
}

SparseMatrix SparseMatrix::getFormat(){
    cout << "look! a wild format!" << endl;
    return SparseMatrix(this->_ncols,this->_nrows);
}



};