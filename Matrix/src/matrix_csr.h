#pragma once

#include "def.h"

#include <iostream>

namespace math
{
  class MatrixCSR
  {
  private:
    size_t _width;
    size_t _height;
    size_t _elementCount;
    size_t _reservedSize;
    size_t _reservedHeight;

    double* _data;
    size_t* _jptr;
    size_t* _iptr;
  public:
    MatrixCSR() noexcept;
    MatrixCSR(const MatrixCSR& matr) noexcept;
    MatrixCSR(MatrixCSR&& matr) noexcept;
    MatrixCSR(const Matrix& matr) noexcept;
    MatrixCSR(size_t width, size_t height) noexcept;
    MatrixCSR(size_t size) noexcept;
    ~MatrixCSR() noexcept;

    MatrixCSR& operator=(const MatrixCSR& matr) noexcept;
    MatrixCSR& operator=(MatrixCSR&& matr) noexcept;
    MatrixCSR& operator=(const Matrix& matr) noexcept;

    MatrixCSR operator*(const MatrixCSR& matr) const noexcept;
    MatrixCSR& operator*=(const MatrixCSR& matr) noexcept;

    operator Matrix() const noexcept;
    operator const Matrix() const noexcept;

    double getElement(size_t y, size_t x) const noexcept;
    double getElement(size_t num) const noexcept;
    void setElement(size_t y, size_t x, double value) noexcept;
    void setElement(size_t num, double value) noexcept;
    double operator[](size_t num) const noexcept;

    void reserve(size_t capacity) noexcept;
    void shrink_to_fit() noexcept;

    size_t size() const noexcept;
    size_t capacity() const noexcept;
    size_t getWidth() const noexcept;
    size_t getHeight() const noexcept;

    void lu(MatrixCSR& l, MatrixCSR& u) const noexcept;
    static MatrixCSR solveSystem(const MatrixCSR& system, const MatrixCSR& answers) noexcept;
    MatrixCSR inverseLU() const noexcept;
    double determinantLU() const noexcept;
    double condLU() const noexcept;

    MatrixCSR& transpose() noexcept;
    MatrixCSR transposing() const noexcept;

    static MatrixCSR identity(size_t size) noexcept;
    static MatrixCSR hilbert(size_t size) noexcept;

    bool isSymmetrical() const noexcept;
    size_t rotatedJacobi(MatrixCSR& coefficients, MatrixCSR& solution, double precision) const noexcept;
    size_t rotateJacobi(MatrixCSR& solution, double precision) noexcept;

  private:
    static MatrixCSR solveSystem(const MatrixCSR& l, const MatrixCSR& u, const MatrixCSR& answers) noexcept;
    double determinantLU(const MatrixCSR& u) const noexcept;
    bool binarySearch(size_t l, size_t r, size_t key, size_t& ind) const noexcept;
  };

  std::ostream& operator<<(std::ostream& out, const MatrixCSR& matrix) noexcept;
}
