#include "matrix_utils.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>

// These are Matrix Constructors
Matrix::Matrix()
    : rows(0), cols(0), data() {}

Matrix::Matrix(int r, int c)
    : rows(r), cols(c), data(r * c, 0.0) {}

Matrix::Matrix(int r, int c, double fill)
    : rows(r), cols(c), data(r * c, fill) {}

// Operator overloading for element access
double &Matrix::operator()(int r, int c)
{
#ifdef DEBUG
    if (r < 0 || r >= rows || c < 0 || c >= cols)
        throw std::out_of_range("Matrix index out of range");
#endif
    return data[r * cols + c];
}

const double &Matrix::operator()(int r, int c) const
{
#ifdef DEBUG
    if (r < 0 || r >= rows || c < 0 || c >= cols)
        throw std::out_of_range("Matrix index out of range");
#endif
    return data[r * cols + c];
}

// Helper Functions
// returns an identity matrix of size nxn
Matrix Matrix::identity(int n)
{
    Matrix I(n, n);
    for (int i = 0; i < n; ++i)
        I(i, i) = 1.0;
    return I;
}

Matrix Matrix::zeros(int r, int c)
{
    return Matrix(r, c); // zero matrix of size rxc
}

// Arithmetic Operations
// Like add, sub, mult, compound operations

Matrix Matrix::operator+(const Matrix &rhs) const
{
    if (rows != rhs.rows || cols != rhs.cols)
        throw std::invalid_argument("Matrix::operator+ : dimension mismatch");
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i)
        result.data[i] = data[i] + rhs.data[i];
    return result;
}

Matrix Matrix::operator-(const Matrix &rhs) const
{
    if (rows != rhs.rows || cols != rhs.cols)
        throw std::invalid_argument("Matrix::operator- : dimension mismatch");
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i)
        result.data[i] = data[i] - rhs.data[i];
    return result;
}

Matrix Matrix::operator*(const Matrix &rhs) const
{
    if (cols != rhs.rows)
        throw std::invalid_argument("Matrix::operator* : dimension mismatch (" + std::to_string(rows) + "x" + std::to_string(cols) + ") * (" + std::to_string(rhs.rows) + "x" + std::to_string(rhs.cols) + ")");

    Matrix result(rows, rhs.cols);

    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < cols; ++k)
        {
            double a_ik = data[i * cols + k];
            for (int j = 0; j < rhs.cols; ++j)
                result.data[i * rhs.cols + j] += a_ik * rhs.data[k * rhs.cols + j];
        }
    return result;
}

Matrix Matrix::operator*(double scalar) const
{
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i)
        result.data[i] = data[i] * scalar;
    return result;
}

Matrix operator*(double scalar, const Matrix &m)
{
    return m * scalar;
}

Matrix &Matrix::operator+=(const Matrix &rhs)
{
    if (rows != rhs.rows || cols != rhs.cols)
        throw std::invalid_argument("Matrix::operator+= : dimension mismatch");
    for (int i = 0; i < rows * cols; ++i)
        data[i] += rhs.data[i];
    return *this;
}

// Transpose of a Matrix
Matrix Matrix::transpose() const
{
    Matrix result(cols, rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result(j, i) = (*this)(i, j);
    return result;
}

// Cholesky decomposition
// Here, we will compute the lower triangular matrix
// It is faster than LU decomposition for SPD matrices and more stable numerically

Matrix Matrix::cholesky() const
{
    if (rows != cols)
        throw std::invalid_argument("cholesky: matrix must be square");

    int n = rows;
    Matrix L(n, n); 

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double sum = 0.0;

            if (j == i)
            {
                for (int k = 0; k < j; ++k)
                    sum += L(j, k) * L(j, k);

                double diag = (*this)(j, j) - sum;
                if (diag <= 0.0)
                    throw std::runtime_error(
                        "cholesky: matrix is not positive definite "
                        "(non-positive pivot at index " +
                        std::to_string(j) + ")");
                L(j, j) = std::sqrt(diag);
            }
            else
            {
                for (int k = 0; k < j; ++k)
                    sum += L(i, k) * L(j, k);
                L(i, j) = ((*this)(i, j) - sum) / L(j, j);
            }
        }
    }
    return L;
}

// Cholesky solve which is A * X = B
// This is how we compute the Kalman gain K = P * Hᵀ * S⁻¹:
// Rearrange to  S * Kᵀ = H * P  then solve for Kᵀ, and then transpose.

Matrix Matrix::cholesky_solve(const Matrix &B) const
{
    if (rows != cols)
        throw std::invalid_argument("cholesky_solve: matrix A must be square");
    if (rows != B.rows)
        throw std::invalid_argument("cholesky_solve: A and B row mismatch");

    int n = rows;
    int m = B.cols;

    Matrix L = this->cholesky();
    Matrix Y(n, m); 
    Matrix X(n, m); 

    for (int col = 0; col < m; ++col)
    {
        for (int i = 0; i < n; ++i)
        {
            double sum = B(i, col);
            for (int k = 0; k < i; ++k)
                sum -= L(i, k) * Y(k, col);
            Y(i, col) = sum / L(i, i);
        }
    }

    for (int col = 0; col < m; ++col)
    {
        for (int i = n - 1; i >= 0; --i)
        {
            double sum = Y(i, col);
            for (int k = i + 1; k < n; ++k)
                sum -= L(k, i) * X(k, col);
            X(i, col) = sum / L(i, i);
        }
    }

    return X;
}

// Block operations
// Useful for block matrices

void Matrix::set_block(int row_off, int col_off, const Matrix &block)
{
    if (row_off + block.rows > rows || col_off + block.cols > cols)
        throw std::out_of_range("set_block: block extends outside matrix bounds");
    for (int i = 0; i < block.rows; ++i)
        for (int j = 0; j < block.cols; ++j)
            (*this)(row_off + i, col_off + j) = block(i, j);
}

Matrix Matrix::get_block(int row_off, int col_off, int r, int c) const
{
    if (row_off + r > rows || col_off + c > cols)
        throw std::out_of_range("get_block: requested block outside matrix bounds");
    Matrix result(r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            result(i, j) = (*this)(row_off + i, col_off + j);
    return result;
}

// Utilities
void Matrix::print(const std::string &label, int max_print) const
{
    if (!label.empty())
        std::cout << "=== " << label << " ===\n";
    std::cout << "Size: " << rows << " x " << cols << "\n";

    int print_rows = std::min(rows, max_print);
    int print_cols = std::min(cols, max_print);

    for (int i = 0; i < print_rows; ++i)
    {
        std::cout << "  [ ";
        for (int j = 0; j < print_cols; ++j)
            std::cout << std::setw(12) << std::setprecision(5)
                      << std::fixed << (*this)(i, j) << " ";
        if (print_cols < cols)
            std::cout << "... ";
        std::cout << "]\n";
    }
    if (print_rows < rows)
        std::cout << "  ... (" << rows - print_rows << " more rows)\n";
    std::cout << "\n";
}

double Matrix::frobenius_norm() const
{
    double sum = 0.0;
    for (double v : data)
        sum += v * v;
    return std::sqrt(sum);
}

bool Matrix::is_symmetric(double tol) const
{
    if (rows != cols)
        return false;
    for (int i = 0; i < rows; ++i)
        for (int j = i + 1; j < cols; ++j)
            if (std::abs((*this)(i, j) - (*this)(j, i)) > tol)
                return false;
    return true;
}

bool Matrix::has_positive_diagonal() const
{
    int n = std::min(rows, cols);
    for (int i = 0; i < n; ++i)
        if ((*this)(i, i) <= 0.0)
            return false;
    return true;
}

void Matrix::symmetrize()
{
    if (rows != cols)
        throw std::invalid_argument("symmetrize: matrix must be square");
    for (int i = 0; i < rows; ++i)
        for (int j = i + 1; j < cols; ++j)
        {
            double avg = 0.5 * ((*this)(i, j) + (*this)(j, i));
            (*this)(i, j) = avg;
            (*this)(j, i) = avg;
        }
}

// Free helper functions
Matrix block_diagonal(const std::vector<Matrix> &blocks)
{
    if (blocks.empty())
        throw std::invalid_argument("block_diagonal: empty block list");

    int total_rows = 0, total_cols = 0;
    for (const auto &b : blocks)
    {
        total_rows += b.rows;
        total_cols += b.cols;
    }

    Matrix result(total_rows, total_cols);

    int row_off = 0, col_off = 0;
    for (const auto &b : blocks)
    {
        result.set_block(row_off, col_off, b);
        row_off += b.rows;
        col_off += b.cols;
    }
    return result;
}

Matrix outer_product(const Matrix &u, const Matrix &v)
{
    if (u.cols != 1 || v.cols != 1)
        throw std::invalid_argument("outer_product: inputs must be column vectors (n×1)");
    Matrix result(u.rows, v.rows);
    for (int i = 0; i < u.rows; ++i)
        for (int j = 0; j < v.rows; ++j)
            result(i, j) = u(i, 0) * v(j, 0);
    return result;
}
