#include<iostream>
#include<vector>
#include<cassert>
#include<cstring>
#include "matrix.h"
using namespace std;

// Initializes matrix without prescribing values
template <class T>
Matrix<T>::Matrix(int nrows, int ncols){
    this->nrows = nrows;
    this->ncols = ncols;
    alloc();
    initialized = true;
}

// Initializes matrix by prescribing `initvalue` to all elements
template <class T>
Matrix<T>::Matrix(int nrows, int ncols, T initvalue){
    this->nrows = nrows;
    this->ncols = ncols;
    alloc();
    initialized = true;
    allfill_(initvalue);
}

// Copy constructor: 
template <class T>
Matrix<T>::Matrix(const Matrix<T>& A){
    this->nrows = A.nrows;
    this->ncols = A.ncols;
    alloc();
    initialized = true;
    for (unsigned int i=0;i<nrows; i++)
        for (unsigned int j=0; j<ncols; j++)
            M[i][j] = A(i, j);
}

template <class T>
Matrix<T>::~Matrix<T>()
{
    allfree();
}

// Fills all elements of matrix by `value`. 
// Assumes that `alloc` has already been called
template <class T>
void Matrix<T>::allfill_(T value){
    assert (this->initialized);
    for (unsigned int i=0; i<this->nrows; i++)
        for (unsigned int j=0; j<this->ncols; j++)
            this->M[i][j] = value;
}

void Matrix<T>::print() const{
    for (unsigned int i=0; i<nrows; i++){
        for (unsigned int j=0; j<ncols; j++)
            cout << this->M[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

// Reference of matrix element x,y
template <class T>
inline T& Matrix<T>::operator()(unsigned int x, unsigned int y) { 
    assert (x <= nrows && y <= ncols);
    return M[x][y]; 
}

template <class T>
inline T Matrix<T>::operator()(unsigned int x, unsigned int y) const { 
    assert (x <= nrows && y <= ncols);
    return M[x][y]; 
}

// Element-wise addition
template <class T>
Matrix& Matrix<T>::operator+=(const Matrix& A)
{
    assert (nrows == A.nrows && ncols == A.ncols);
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] += A.M[i][j];
    return *this;
}

// Scalar multiplication
template <class T>
Matrix& Matrix<T>::operator*=(T c)
{
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] *= c;
    return *this;
}

template <class T>
Matrix& Matrix<T>::operator-=(const Matrix& A)
{
    assert (nrows == A.nrows && ncols == A.ncols);
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] -= A.M[i][j];
    return *this;
}

template <class T>
Matrix& Matrix<T>::operator=(const Matrix& A){
    if (this == &A)
        return *this;
    allfree();
    nrows = A.nrows;
    ncols = A.ncols;
    alloc();
    initialized = true;
    for (unsigned int i=0;i<nrows; i++)
        for (unsigned int j=0; j<ncols; j++)
            M[i][j] = A(i, j);
    return *this;
}

// Provides copy of matrix slice A[ri:rj, ci:cj] supporting negative indices
template <class T>
Matrix Matrix<T>::slice(int ri, int rj, int ci, int cj) const{
    rj = rj < 0 ? nrows + rj : rj;
    cj = cj < 0 ? ncols + cj : cj;
    ri = ri < 0 ? nrows + ri : ri;
    ci = ci < 0 ? ncols + ci : ci;
    assert (ri >= 0 && ri < rj && (unsigned int) rj <= nrows && ci >= 0 && ci < cj && (unsigned int) cj <= ncols);
    Matrix<T> B(rj-ri, cj-ci);
    for (unsigned int i=0; i<B.nrows; ++i)
        for (unsigned int j=0; j<B.ncols; ++j)
            B(i, j) = M[i+ri][j+ci];
    return B;
}

// Appends B to current matrix by axis
template <class T>
Matrix Matrix<T>::append(const Matrix<T>& B, unsigned int axis) const {
    assert (axis == 0 || axis == 1);
    unsigned int nr = nrows + (1-axis) * B.nrows, nc = ncols + axis * B.ncols; 
    Matrix<T> C(nr, nc);

    for (unsigned int i=0; i<nr; i++)
        for (unsigned int j=0; j<nc; j++){
            if (i >= nrows){
                C(i, j) = B(i-nrows, j);
                continue;
            }
            if (j >= ncols){
                C(i, j) = B(i, j-ncols);
                continue;
            }
            C(i, j) = M[i][j];
        }
    return C;
}

template <class T>
void Matrix<T>::alloc(){
    M = new T*[nrows];
    for (unsigned int i=0; i<nrows; i++)
        M[i] = new T[ncols];
}

template <class T>
void Matrix<T>::allfree(){
    for (unsigned int i = 0; i < this->nrows; ++i) {
        delete[] M[i];
    }
    delete[] M;
}