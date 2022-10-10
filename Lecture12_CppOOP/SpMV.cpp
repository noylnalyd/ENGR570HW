#include <iostream>
#include "SpMV.hpp"


int main(int argc, char const *argv[])
{
    
    // allocate spmv
    SpMV::SparseMatrix* ptr_A = new SpMV::SparseMatrix(2,5);

    // Scoping unit
    {
        SpMV::SparseMatrix A = SpMV::SparseMatrix(3,4);
    }

    delete(ptr_A);
    ptr_A = NULL;
    return 0;
}
