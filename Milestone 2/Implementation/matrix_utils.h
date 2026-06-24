#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H
// The purpose is to provide a heap-allocated Matrix class for all Kalman Filter operations.
// All matrices are allocated on the HEAP using std::vector<double>.
// The reason we are doing this is that the full-body state space requires matrices as large as 276×276
// and so Heap allocation also allows matrix sizes to be determined at runtime,
// making the code flexible for any number of joints.
#include <vector>
#include <string>
#include <stdexcept>

// Matrix class
class Matrix {
public:
    int rows;
    int cols;
    std::vector<double> data; 

    // constructors
    Matrix();

    Matrix(int rows, int cols);

    Matrix(int rows, int cols, double fill_value);

    // element access
    double&       operator()(int r, int c);
    const double& operator()(int r, int c) const;

    // n×n identity matrix
    static Matrix identity(int n);

    // r×c zero matrix
    static Matrix zeros(int r, int c);

    // Arithmetic operators

    Matrix operator+(const Matrix& rhs) const;
    Matrix operator-(const Matrix& rhs) const;
    Matrix operator*(const Matrix& rhs) const;

    // Scalar multiplication
    Matrix operator*(double scalar) const;
    friend Matrix operator*(double scalar, const Matrix& m);

    // In-place addition
    Matrix& operator+=(const Matrix& rhs);

    // Transpose 
    Matrix transpose() const;

    // Cholesky decomposition which would return a lower-triangular L such that A = L * Lᵀ
    Matrix cholesky() const;

    // Solve  A * X = B 
    Matrix cholesky_solve(const Matrix& B) const;

    // Block operations
    // Copy block into this matrix starting at (row_off, col_off)
    void set_block(int row_off, int col_off, const Matrix& block);

    // Extract sub-matrix of size (r × c) starting at (row_off, col_off)
    Matrix get_block(int row_off, int col_off, int r, int c) const;

    // Utilities
    void print(const std::string& label = "", int max_print = 8) const;

    // Frobenius norm which would be useful in debugging
    double frobenius_norm() const;

    // Check if matrix is approx symmetric
    bool is_symmetric(double tol = 1e-9) const;

    // Returns true if all diagonal elements are positive. This is a SPD check
    bool has_positive_diagonal() const;

    // Symmetrize:  A = (A + Aᵀ) / 2 
    void symmetrize();
};

// Other Helper functions
// Build a block-diagonal matrix from a vector of same-size blocks.
// All blocks must have identical dimensions.
Matrix block_diagonal(const std::vector<Matrix>& blocks);

Matrix outer_product(const Matrix& u, const Matrix& v);

#endif 
