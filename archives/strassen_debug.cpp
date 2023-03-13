#include<iostream>
#include<vector>
#include<cassert>
#include<cstring>
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
        Matrix(int nrows, int ncols){
            this->nrows = nrows;
            this->ncols = ncols;
            alloc();
            initialized = true;
        }

        // Initializes matrix by prescribing `initvalue` to all elements
        Matrix(int nrows, int ncols, T initvalue){
            this->nrows = nrows;
            this->ncols = ncols;
            alloc();
            initialized = true;
            allfill_(initvalue);
        }

        // Initializes matrix by prescribing `initvalue` to all elements
        Matrix(const Matrix<T>& A){
            cout << "Copy constructor called! on "<< endl;
            A.print();
            this->nrows = A.nrows;
            this->ncols = A.ncols;
            alloc();
            initialized = true;
            for (unsigned int i=0;i<nrows; i++)
                for (unsigned int j=0; j<ncols; j++)
                    M[i][j] = A(i, j);
            cout << "Copy ended" << endl;
        }

        Matrix(Matrix<T>&& A){
            cout << "Move constructor called! on "<< endl;
            A.print();
            this->nrows = A.nrows;
            this->ncols = A.ncols;
            alloc();
            initialized = true;
            for (unsigned int i=0;i<nrows; i++)
                for (unsigned int j=0; j<ncols; j++)
                    M[i][j] = A(i, j);
            this->print();
            cout << "Copy ended" << endl;
        }

        ~Matrix()
        {
            allfree();
        }

        // Fills all elements of matrix by `value`. 
        // Assumes that `alloc` has already been called
        void allfill_(T value){
            assert (this->initialized);
            for (unsigned int i=0; i<this->nrows; i++)
                for (unsigned int j=0; j<this->ncols; j++)
                    this->M[i][j] = value;
        }

        void print() const{
            for (unsigned int i=0; i<nrows; i++){
                for (unsigned int j=0; j<ncols; j++)
                    cout << this->M[i][j] << " ";
                cout << endl;
            }
            cout << endl;
        }

        // Reference of matrix element x,y
        inline T& operator()(unsigned int x, unsigned int y) { 
            assert (x <= nrows && y <= ncols);
            return M[x][y]; 
        }

        inline T operator()(unsigned int x, unsigned int y) const { 
            assert (x <= nrows && y <= ncols);
            return M[x][y]; 
        }

        // Element-wise addition
        Matrix& operator+=(const Matrix& A)
        {
            assert (nrows == A.nrows && ncols == A.ncols);
            for (unsigned int i = 0; i < nrows; ++i) 
                for (unsigned int j = 0; j < ncols; ++j) 
                    M[i][j] += A.M[i][j];
            return *this;
        }

        // Scalar multiplication
        Matrix& operator*=(T c)
        {
            for (unsigned int i = 0; i < nrows; ++i) 
                for (unsigned int j = 0; j < ncols; ++j) 
                    M[i][j] *= c;
            return *this;
        }

        Matrix& operator-=(const Matrix& A)
        {
            assert (nrows == A.nrows && ncols == A.ncols);
            for (unsigned int i = 0; i < nrows; ++i) 
                for (unsigned int j = 0; j < ncols; ++j) 
                    M[i][j] -= A.M[i][j];
            return *this;
        }

        Matrix& operator=(const Matrix& A){
            cout << "Copy Called!" << endl;
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
        Matrix slice(int ri, int rj, int ci, int cj) const{
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
        Matrix append(const Matrix<T>& B, unsigned int axis) const {
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
    private:
        void alloc(){
            M = new T*[nrows];
            for (unsigned int i=0; i<nrows; i++)
                M[i] = new T[ncols];
        }

        void allfree(){
            for (unsigned int i = 0; i < this->nrows; ++i) {
                delete[] M[i];
            }
            delete[] M;
        }
};

template<class T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B){
    Matrix<T> C(A);
    return (C += B);
}

template<class T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B){
    Matrix<T> C(A.nrows, A.ncols);
    C += A;
    C -= B;
    return C;
}

template<class T>
Matrix<T> operator*(const Matrix<T>& A, T c){
    Matrix<T> C(A.nrows, A.ncols);
    C += A;
    C *= c;
    return C;
}

// Naive algorithm for matrix multiplication
template<class T>
Matrix<T> simplematmul(Matrix<T>& A, Matrix<T>& B){
    assert (A.ncols == B.nrows);
    Matrix<T> C(A.nrows, A.ncols, (T)(0));
    for (unsigned int i=0; i<A.nrows; i++)
        for (unsigned int k=0; k<A.ncols; k++)
            for (unsigned int j=0; j<B.ncols; j++)
                C(i, j) += A(i, k) * B(k, j);
    return C;
}

template<class T>
// Pads matrix to square matrix with even dimensions
// Always returns copy
Matrix<T> pad(const Matrix<T>& A){
    assert (A.nrows == A.ncols);
    if (A.nrows % 2 == 0)
        return A.slice(0, A.nrows, 0, A.ncols);
    return A.append(*(new Matrix<T>(1, A.ncols, (T)0)), 0).append(*(new Matrix<T>(A.nrows+1, 1, (T)0)), 1);
}

template<class T>
Matrix<T> dot(const Matrix<T>& u, const Matrix<T>& v){
    assert (u.ncols == v.nrows && u.nrows == 1 && v.ncols == 1);
    T a = (T) 0;
    for (unsigned int i=0;i<u.nrows;i++)
        a += u(0, i) * v(i, 0);
    return *(new Matrix<T>(1, 1, a));
}

template<class T>
Matrix<T> strassen(const Matrix<T>& X, const Matrix<T>& Y){
    assert (X.ncols == Y.nrows);
    X.print();
    Y.print();
    if (X.nrows == 1 && Y.ncols == 1)
        return dot(X, Y);
    // We only implement for square matrices
    if (X.ncols % 2 != 0 || X.nrows % 2 != 0 || Y.ncols % 2 != 0){
        return strassen(pad(X), pad(Y)).slice(0, -1, 0, -1);
    }
    Matrix<T> A = X.slice(0, X.nrows / 2, 0, X.ncols / 2), B = X.slice(0, X.nrows / 2, X.ncols / 2, X.ncols), 
              C = X.slice(X.nrows / 2, X.nrows, 0, X.ncols / 2), D = X.slice(X.nrows / 2, X.nrows, X.ncols / 2, X.ncols),
              E = Y.slice(0, Y.nrows / 2, 0, Y.ncols / 2), F = Y.slice(0, Y.nrows / 2, Y.ncols / 2, Y.ncols), 
              G = Y.slice(Y.nrows / 2, Y.nrows, 0, Y.ncols / 2), H = Y.slice(Y.nrows / 2, Y.nrows, Y.ncols / 2, Y.ncols);
    A.print();
    (F-H).print();
    auto x1 = F - H;
    auto P1 = strassen(A, x1), P2 = strassen(A+B, H), P3 = strassen(C+D, E), P4 = strassen(D, G-E), 
            P5 = strassen(A+D, E+H), P6 = strassen(B-D, G+H), P7 = strassen(C-A, E+F);
    return ((P4 + P5 + P6 - P2).append(P1+P2, 1)).append((P3+P4).append(P1 - P3 + P5 + P7, 1), 0);
}

int main(){
    unsigned int D = 2;
    Matrix<int> A(D, D, -2), B(D, D, 1);
    auto C = A + B;
}