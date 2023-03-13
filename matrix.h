#ifndef _MATRIX_H_
#define _MATRIX_H_

#include<iostream>
#include<cassert>
using namespace std;

template <class T>
class Matrix{
    public:
        T** M;
        // Effective matrix is M[rs:rs+sz, cs:cs+sz]
        unsigned int nrows;
        unsigned int ncols; 
        bool initialized = false; 

        // Initializes matrix without prescribing values
        Matrix(int nrows, int ncols);

        // Initializes matrix by prescribing `initvalue` to all elements
        Matrix(int nrows, int ncols, T initvalue);

        // Copy constructor: 
        Matrix(const Matrix<T>& A);

        ~Matrix();

        // Fills all elements of matrix by `value`. 
        // Assumes that `alloc` has already been called
        void allfill_(T value);

        void print() const;

        // Reference of matrix element x,y
        inline T& operator()(unsigned int x, unsigned int y);

        // Element-wise addition
        Matrix& operator+=(const Matrix& A);

        // Scalar multiplication
        Matrix& operator*=(T c);

        Matrix& operator-=(const Matrix& A);

        Matrix& operator=(const Matrix& A);

        // Provides copy of matrix slice A[ri:rj, ci:cj] supporting negative indices
        Matrix slice(int ri, int rj, int ci, int cj) const;

        // Appends B to current matrix by axis
        Matrix append(const Matrix<T>& B, unsigned int axis) const;

        Matrix simplematmul(const Matrix& B);
        Matrix strassen(const Matrix&, unsigned int);
        Matrix pad() const;
    private:
        void alloc();
        void allfree();
};

template<class T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B){
    Matrix<T> C(A);
    return (C += B);
}
template<class T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B){
    Matrix<T> C(A);
    return (C -= B);
}

template<class T>
Matrix<T> operator*(const Matrix<T>& A, T c){
    Matrix<T> C(A);
    return (C *= c);
}
template<class T>
Matrix<T> dot(const Matrix<T>& u, const Matrix<T>& v){
    assert (u.ncols == v.nrows && u.nrows == 1 && v.ncols == 1);
    T a = (T) 0;
    for (unsigned int i=0;i<u.nrows;i++)
        a += u.M[0][i] * v.M[i][0];
    return *(new Matrix<T>(1, 1, a));
}


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
            M[i][j] = A.M[i][j];
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

template <class T>
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

// Element-wise addition
template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& A)
{
    assert (nrows == A.nrows && ncols == A.ncols);
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] += A.M[i][j];
    return *this;
}

// Scalar multiplication
template <class T>
Matrix<T>& Matrix<T>::operator*=(T c)
{
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] *= c;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& A)
{
    assert (nrows == A.nrows && ncols == A.ncols);
    for (unsigned int i = 0; i < nrows; ++i) 
        for (unsigned int j = 0; j < ncols; ++j) 
            M[i][j] -= A.M[i][j];
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& A){
    if (this == &A)
        return *this;
    allfree();
    nrows = A.nrows;
    ncols = A.ncols;
    alloc();
    initialized = true;
    for (unsigned int i=0;i<nrows; i++)
        for (unsigned int j=0; j<ncols; j++)
            M[i][j] = A.M[i][j];
    return *this;
}

// Provides copy of matrix slice A[ri:rj, ci:cj] supporting negative indices
template <class T>
Matrix<T> Matrix<T>::slice(int ri, int rj, int ci, int cj) const{
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
Matrix<T> Matrix<T>::append(const Matrix<T>& B, unsigned int axis) const {
    assert (axis == 0 || axis == 1);
    unsigned int nr = nrows + (1-axis) * B.nrows, nc = ncols + axis * B.ncols; 
    Matrix<T> C(nr, nc);

    for (unsigned int i=0; i<nr; i++)
        for (unsigned int j=0; j<nc; j++){
            if (i >= nrows){
                C(i, j) = B.M[i-nrows][j];
                continue;
            }
            if (j >= ncols){
                C(i, j) = B.M[i][j-ncols];
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

template <class T>
Matrix<T> Matrix<T>::simplematmul(const Matrix& B){
    assert (ncols == nrows);
    Matrix<T> C(nrows, ncols, (T)(0));
    for (unsigned int i=0; i<nrows; i++)
        for (unsigned int k=0; k<ncols; k++)
            for (unsigned int j=0; j<B.ncols; j++)
                C(i, j) += M[i][k] * B.M[k][j];
    return C;
}

template <class T>
Matrix<T> Matrix<T>::pad() const{
    assert (nrows == ncols);
    if (nrows % 2 == 0){
        auto t = *this;
        return t;
    }
    return this->append(*(new Matrix<T>(1, ncols, (T)0)), 0).append(*(new Matrix<T>(nrows+1, 1, (T)0)), 1);
}

template <class T>
Matrix<T> Matrix<T>::strassen(const Matrix<T>& Y, unsigned int threshold){
    assert (ncols == Y.nrows);
    if (nrows == 1 && Y.ncols == 1)
        return dot(*this, Y);
    // We only implement for square matrices
    if (ncols % 2 != 0 || nrows % 2 != 0 || Y.ncols % 2 != 0){
        return (this->pad()).strassen(Y.pad(), threshold).slice(0, -1, 0, -1);
    }
    Matrix<T> A = this->slice(0, nrows / 2, 0, ncols / 2), B = this->slice(0, nrows / 2, ncols / 2, ncols), 
              C = this->slice(nrows / 2, nrows, 0, ncols / 2), D = this->slice(nrows / 2, nrows, ncols / 2, ncols),
              E = Y.slice(0, Y.nrows / 2, 0, Y.ncols / 2), F = Y.slice(0, Y.nrows / 2, Y.ncols / 2, Y.ncols), 
              G = Y.slice(Y.nrows / 2, Y.nrows, 0, Y.ncols / 2), H = Y.slice(Y.nrows / 2, Y.nrows, Y.ncols / 2, Y.ncols);
    auto x1 = F - H;
    auto P1 = A.strassen(x1, threshold), P2 = (A+B).strassen(H, threshold), P3 = (C+D).strassen(E, threshold), P4 = D.strassen(G-E, threshold), 
            P5 = (A+D).strassen(E+H, threshold), P6 = (B-D).strassen(G+H, threshold), P7 = (C-A).strassen(E+F, threshold);
    return ((P4 + P5 + P6 - P2).append(P1+P2, 1)).append((P3+P4).append(P1 - P3 + P5 + P7, 1), 0);
}
#endif