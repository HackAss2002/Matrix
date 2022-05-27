#include "matrix_csr.h"
#include "matrix.h"

#include <cstring>
#include <cassert>

namespace math
{

Matrix::Matrix() noexcept : _width(0), _height(0), _reservedSize(0), _data(nullptr)
{

}

Matrix::Matrix(const Matrix& matr) noexcept : _width(matr._width), _height(matr._height),
  _reservedSize(matr._width * matr._height), _data(new double[matr._height * matr._width])
{
  memcpy(_data, matr._data, sizeof(double) * _reservedSize);
}

Matrix::Matrix(Matrix&& matr) noexcept : _width(matr._width), _height(matr._height),
  _reservedSize(matr._reservedSize), _data(matr._data)
{
  matr._width = 0;
  matr._height = 0;
  matr._reservedSize = 0;
  matr._data = nullptr;
}

Matrix::Matrix(const MatrixCSR& matr) noexcept : _width(matr.getWidth()), _height(matr.getHeight()),
  _reservedSize(matr.getWidth()* matr.getHeight()), _data(new double[matr.getWidth() * matr.getHeight()])
{
  for (size_t num = 0; num < _reservedSize; ++num)
    _data[num] = matr[num];
}

Matrix::Matrix(size_t width, size_t height) noexcept : _width(width), _height(height),
  _reservedSize(width * height), _data(new double[width * height])
{
  memset(_data, 0, sizeof(double) * _reservedSize);
}

Matrix::Matrix(size_t size) noexcept : _width(size), _height(size),
  _reservedSize(size * size), _data(new double[size * size])
{
  memset(_data, 0, sizeof(double) * _reservedSize);
}

Matrix::~Matrix() noexcept
{
  delete[] _data;
}

Matrix& Matrix::operator=(const Matrix& matr) noexcept
{
  if (this == &matr)
    return *this;

  _width = matr._width;
  _height = matr._height;

  const size_t elementCount = _width * _height;
  if (elementCount > _reservedSize)
  {
    double* newData = new double[elementCount];
    memcpy(newData, matr._data, sizeof(double) * elementCount);
    delete[] _data;
    _data = newData;
    _reservedSize = elementCount;
  }
  else
    memcpy(_data, matr._data, sizeof(double) * elementCount);

  return *this;
}

Matrix& Matrix::operator=(Matrix&& matr) noexcept
{
  if (this == &matr)
    return *this;

  _width = matr._width;
  _height = matr._height;
  _reservedSize = matr._reservedSize;

  matr._width = 0;
  matr._height = 0;
  matr._reservedSize = 0;

  delete[] _data;
  _data = matr._data;
  matr._data = nullptr;

  return *this;
}

Matrix& Matrix::operator=(const MatrixCSR& matr) noexcept
{
  _width = matr.getWidth();
  _height = matr.getHeight();

  const size_t elementCount = _width * _height;
  if (elementCount > _reservedSize)
  {
    delete[] _data;
    _data = new double[elementCount];
    _reservedSize = elementCount;
  }

  for (size_t num = 0; num < elementCount; ++num)
    _data[num] = matr[num];

  return *this;
}

Matrix Matrix::operator*(const Matrix& matr) const noexcept
{
  assert(_width == matr._height);
  size_t i, j, k;
  double val;
  Matrix res(matr._width, _height);
  for (i = 0; i < res._height; ++i)
  {
    const size_t offset = i * _width;
    for (j = 0; j < res._width; ++j)
    {
      val = 0;
      for (k = 0; k < _width; ++k)
        val += (*this)[offset + k] * matr[k * matr._width + j];
      res.setElement(i, j, val);
    }
  }
  return res;
}

Matrix& Matrix::operator*=(const Matrix& matr) noexcept
{
  *this = *this * matr;
  return *this;
}

Matrix::operator MatrixCSR() const noexcept
{
  const size_t elementCount = _width * _height;
  size_t countNotZero = 0;
  for (size_t num = 0; num < elementCount; ++num)
    countNotZero += !!_data[num];
  MatrixCSR matr;
  matr.reserve(countNotZero);
  for (size_t num = 0; num < elementCount; ++num)
    if (_data[num])
      matr.setElement(num, _data[num]);
  return matr;
}

Matrix::operator const MatrixCSR() const noexcept
{
  const size_t elementCount = _width * _height;
  size_t countNotZero = 0;
  for (size_t num = 0; num < elementCount; ++num)
    countNotZero += !!_data[num];
  MatrixCSR matr;
  matr.reserve(countNotZero);
  for (size_t num = 0; num < elementCount; ++num)
    if (_data[num])
      matr.setElement(num, _data[num]);
  return matr;
}

double Matrix::getElement(size_t y, size_t x) const noexcept
{
  assert(y < _height);
  assert(x < _width);
  return getElement(y * _width + x);
}

double Matrix::getElement(size_t num) const noexcept
{
  assert(num < _width * _height);
  return _data[num];
}

void Matrix::setElement(size_t y, size_t x, double value) noexcept
{
  assert(y < _height);
  assert(x < _width);
  setElement(y * _width + x, value);
}

void Matrix::setElement(size_t num, double value) noexcept
{
  assert(num < _width * _height);
  _data[num] = value;
}

double Matrix::operator[](size_t num) const noexcept
{
  assert(num < _width* _height);
  return _data[num];
}

double& Matrix::operator[](size_t num) noexcept
{
  assert(num < _width* _height);
  return _data[num];
}

Matrix::operator double* () noexcept
{
  return _data;
}

Matrix::operator double const* () const noexcept
{
  return _data;
}

void Matrix::reserve(size_t capacity) noexcept
{
  if (capacity <= _reservedSize)
    return;

  double* newData = new double[capacity];
  memcpy(newData, _data, sizeof(double) * capacity);
  delete[] _data;
  _data = newData;
  _reservedSize = capacity;
}

void Matrix::shrink_to_fit() noexcept
{
  const size_t elementCount = _width * _height;
  if (elementCount == _reservedSize)
    return;

  if (elementCount)
  {
    double* newData = new double[elementCount];
    memcpy(newData, _data, sizeof(double) * elementCount);
    delete[] _data;
    _data = newData;
  }
  else
  {
    delete[] _data;
    _data = nullptr;
  }

  _reservedSize = elementCount;
}

size_t Matrix::size() const noexcept
{
  return _width * _height;
}

size_t Matrix::capacity() const noexcept
{
  return _reservedSize;
}

size_t Matrix::getWidth() const noexcept
{
  return _width;
}

size_t Matrix::getHeight() const noexcept
{
  return _height;
}

void Matrix::lu(Matrix& l, Matrix& u) const noexcept
{
  assert(_width == _height);
  u = *this;
  l = Matrix(_width);

  size_t i, j, k;
  double val;
  for (k = 1; k < _width; ++k)
  {
    for (i = k - 1; i < _width; ++i)
    {
      val = u.getElement(i, i);
      for (j = i; j < _width; ++j)
        l.setElement(j, i, u.getElement(j, i) / val);
    }

    for (i = k; i < _width; ++i)
    {
      val = l.getElement(i, k - 1);
      for (j = k - 1; j < _width; ++j)
        u.setElement(i, j, u.getElement(i, j) - val * u.getElement(k - 1, j));
    }
  }
}

Matrix Matrix::solveSystem(const Matrix& system, const Matrix& answers) noexcept
{
  assert(system._width == system._height);
  assert(system._height == answers._height);
  assert(answers._width == 1);

  Matrix l, u;
  system.lu(l, u);
  return solveSystem(l, u, answers);
}

Matrix Matrix::solveSystem(const Matrix& l, const Matrix& u, const Matrix& answers) noexcept
{
  assert(l._width == l._height);
  assert(l._width == u._width);
  assert(u._width == u._height);
  assert(l._height == answers._height);
  assert(answers._width == 1);

  Matrix tmpRes(1, l._width);
  for (size_t i = 0; i < tmpRes.getHeight(); ++i)
  {
    double sum = 0;
    for (size_t j = 0; j < i; ++j)
      sum += l.getElement(i, j) * tmpRes[j];
    tmpRes.setElement(i, answers[i] - sum);
  }

  Matrix res(1, l._width);
  for (size_t i = res.getHeight(); i > 0; --i)
  {
    double sum = 0;
    for (size_t j = res.getHeight(); j > i; --j)
      sum += u.getElement(i - 1, j - 1) * res[j - 1];
    res.setElement(i - 1, (tmpRes[i - 1] - sum) / u.getElement(i - 1, i - 1));
  }
  return res;
}

Matrix Matrix::inverseLU() const noexcept
{
  assert(_width == _height);
  Matrix inversed(_width, _height);
  Matrix ans(1, _width);
  Matrix l, u;
  lu(l, u);
  size_t i, j;
  for (i = 0; i < _width; ++i)
  {
    ans.setElement(i, 1);
    const Matrix solved = solveSystem(l, u, ans);
    for (j = 0; j < _height; ++j)
      inversed.setElement(j, i, solved[j]);
    ans.setElement(i, 0);
  }
  return inversed;
}

double Matrix::determinantLU() const noexcept
{
  Matrix l, u;
  lu(l, u);
  return determinantLU(u);
}

double Matrix::determinantLU(const Matrix& u) const noexcept
{
  assert(u._width == u._height);
  double res = 1;
  for (size_t i = 0; i < u._width; ++i)
    res *= getElement(i, i);
  return res;
}

double Matrix::condLU() const noexcept
{
  double norm = 0;
  double elem;
  for (size_t i = 0; i < _height; ++i)
  {
    for (size_t j = 0; j < _width; ++j)
    {
      elem = getElement(i, j);
      norm += elem * elem;
    }
  }
  const Matrix inversed = inverseLU();
  double normInv = 0;
  for (size_t i = 0; i < _height; ++i)
  {
    for (size_t j = 0; j < _width; ++j)
    {
      elem = inversed.getElement(i, j);
      normInv += elem * elem;
    }
  }
  return sqrt(norm * normInv);
}

Matrix& Matrix::transpose() noexcept
{
  *this = transposing();
  return *this;
}

Matrix Matrix::transposing() const noexcept
{
  Matrix result(_height, _width);
  size_t i, j;
  for (i = 0; i < result._height; ++i)
    for (j = 0; j < result._width; ++j)
      result.setElement(i, j, getElement(j, i));
  return result;
}

Matrix Matrix::identity(size_t size) noexcept
{
  Matrix matr(size);
  for (size_t i = 0; i < size; ++i)
    matr.setElement(i, i, 1);
  return matr;
}

Matrix Matrix::hilbert(size_t size) noexcept
{
  Matrix matr(size);
  size_t i, j;
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      matr.setElement(i, j, 1 / static_cast<double>(i + j + 1));
  return matr;
}

bool Matrix::isSymmetrical() const noexcept
{
  if (_width != _height)
    return false;

  for (size_t i = 0; i < _width; ++i)
    for (size_t j = 0; j < i; ++j)
      if (getElement(i, j) != getElement(j, i))
        return false;

  return true;
}

size_t Matrix::rotateJacobi(Matrix& solution, double precision) noexcept
{
  assert(isSymmetrical());
  size_t result = 0;
  size_t i, j, k;
  size_t maxI, maxJ;
  double max, fi;
  Matrix matricaPoworota(_width);
  Matrix temp(_width);
  solution = identity(_width);
  double fault = 0.0;
  double value;
  double sqr2 = sqrt(2);
  double precision2 = precision * precision / 2;
  for (i = 0; i < _width; ++i)
    for (j = i + 1; j < _width; ++j)
    {
      value = getElement(i, j);
      fault += value * value;
    }

  while (fault > precision2)
  {
    max = 0.0;
    for (i = 0; i < _width; ++i)
      for (j = i + 1; j < _width; ++j)
      {
        value = fabs(getElement(i, j));
        if (value > max)
        {
          max = value;
          maxI = i;
          maxJ = j;
        }
      }

    matricaPoworota = identity(_width);
    if (getElement(maxI, maxI) == getElement(maxJ, maxJ))
    {
      value = 1 / sqr2;
      matricaPoworota.setElement(maxI, maxI, value);
      matricaPoworota.setElement(maxJ, maxJ, value);
      matricaPoworota.setElement(maxJ, maxI, value);
      matricaPoworota.setElement(maxI, maxJ, -value);
    }
    else
    {
      fi = 0.5 * atan((2.0 * getElement(maxI, maxJ)) /
        (getElement(maxI, maxI) - getElement(maxJ, maxJ)));
      value = cos(fi);
      matricaPoworota.setElement(maxI, maxI, value);
      matricaPoworota.setElement(maxJ, maxJ, value);
      value = sin(fi);
      matricaPoworota.setElement(maxJ, maxI, value);
      matricaPoworota.setElement(maxI, maxJ, -value);
    }

    for (i = 0; i < _width; ++i)
      for (j = 0; j < _width; ++j)
      {
        value = 0;
        for (k = 0; k < _width; ++k)
          value += matricaPoworota.getElement(k, i) * getElement(k, j);
        temp.setElement(i, j, value);
      }

    *this = temp * matricaPoworota;
    fault = 0.0;
    for (i = 0; i < _width; ++i)
      for (j = i + 1; j < _width; ++j)
      {
        value = getElement(i, j);
        fault += value * value;
      }

    solution *= matricaPoworota;
    ++result;
  }
  return result;
}

size_t Matrix::rotatedJacobi(Matrix& coef, Matrix& solution, double precision) const noexcept
{
  coef = *this;
  
  return coef.rotateJacobi(solution, precision);
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) noexcept
{
  out << std::endl;
  for (size_t y = 0; y < matrix.getHeight(); ++y)
  {
    for (size_t x = 0; x < matrix.getWidth(); ++x)
      out << matrix.getElement(y, x) << " ";
    out << std::endl;
  }
  return out;
}

}
