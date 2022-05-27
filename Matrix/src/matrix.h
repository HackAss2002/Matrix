#pragma once

#include "def.h"

#include <iostream>

namespace math
{
  class Matrix
  {
  private:
    size_t _width;
    size_t _height;
    double* _data;
    size_t _reservedSize;
  public:
    Matrix() noexcept;
    Matrix(const Matrix& matr) noexcept;
    Matrix(Matrix&& matr) noexcept;
    Matrix(const MatrixCSR& matr) noexcept;
    Matrix(size_t width, size_t height) noexcept;
    Matrix(size_t size) noexcept;
    ~Matrix() noexcept;

    Matrix& operator=(const Matrix& matr) noexcept;
    Matrix& operator=(Matrix&& matr) noexcept;
    Matrix& operator=(const MatrixCSR& matr) noexcept;

    Matrix operator*(const Matrix& matr) const noexcept;
    Matrix& operator*=(const Matrix& matr) noexcept;

    operator MatrixCSR() const noexcept;
    operator const MatrixCSR() const noexcept;

    double getElement(size_t y, size_t x) const noexcept;
    double getElement(size_t num) const noexcept;
    void setElement(size_t y, size_t x, double value) noexcept;
    void setElement(size_t num, double value) noexcept;
    double operator[](size_t num) const noexcept;
    double& operator[](size_t num) noexcept;

    explicit operator double* () noexcept;
    explicit operator double const* () const noexcept;

    void reserve(size_t capacity) noexcept;
    void shrink_to_fit() noexcept;

    size_t size() const noexcept;
    size_t capacity() const noexcept;
    size_t getWidth() const noexcept;
    size_t getHeight() const noexcept;

    void lu(Matrix& l, Matrix& u) const noexcept;
    static Matrix solveSystem(const Matrix& system, const Matrix& answers) noexcept;
    Matrix inverseLU() const noexcept;
    double determinantLU() const noexcept;
    double condLU() const noexcept;

    Matrix& transpose() noexcept;
    Matrix transposing() const noexcept;

    static Matrix identity(size_t size) noexcept;
    static Matrix hilbert(size_t size) noexcept;

    bool isSymmetrical() const noexcept;
    size_t rotatedJacobi(Matrix& coefficients, Matrix& solution, double precision) const noexcept;
    size_t rotateJacobi(Matrix& solution, double precision) noexcept;

  private:
    static Matrix solveSystem(const Matrix& l, const Matrix& u, const Matrix& answers) noexcept;
    double determinantLU(const Matrix& u) const noexcept;
  };

  std::ostream& operator<<(std::ostream& out, const Matrix& matrix) noexcept;
}
