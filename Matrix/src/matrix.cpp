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

Matrix Matrix::rotatedJacobi(Matrix& v, Matrix& d) const noexcept
{
  assert(_width == _height);
  size_t j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, aipq, dip, diq;
  Matrix b(1, _width), z(1, _width);
  Matrix a = *this;
  v = identity(_width);
  d = Matrix(1, _width);
  const size_t MAXSWEEP = 10000;
  for (ip = 0; ip < _width; ++ip)
    b.setElement(ip, a.getElement(ip, ip));
  size_t nrot = 0;
  auto rotate = [&](Matrix& a, size_t i, size_t j, size_t k, size_t l)
  {
    g = a.getElement(i, j);
    h = a.getElement(k, l);
    a.setElement(i, j, g - s * (h + g * tau));
    a.setElement(k, l, h + s * (g - h * tau));
  };
  for (i = 0; i < MAXSWEEP; ++i)
  {
    for (sm = 0., ip = 0; ip < _width; ++ip)
      for (iq = ip + 1; iq < _width; ++iq)
        sm += fabs(a.getElement(ip, iq));
    if (!sm)
      return a;
    tresh = (i < 3 ? 0.2 * sm / (_width * _width) : 0.);
    for (ip = 0; ip < _width - 1; ++ip)
      for (iq = ip + 1; iq < _width; ++iq)
      {
        aipq = a.getElement(ip, iq);
        g = 100. * fabs(aipq);
        dip = d.getElement(ip);
        diq = d.getElement(iq);
        if (i > 3 && fabs(dip + g) == fabs(dip) && fabs(diq + g) == fabs(diq))
          a.setElement(ip, iq, 0.);
        else if (fabs(aipq) > tresh)
        {
          h = dip - diq;
          if ((fabs(h) + g) == fabs(h))
            t = aipq / h;
          else
          {
            theta = 0.5 * h / aipq;
            t = 1. / (fabs(theta) + sqrt(1. + theta * theta));
            if (theta < 0.)
              t = -t;
          }
          c = 1. / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1. + c);
          h = t * aipq;

          z.setElement(ip, z.getElement(ip) - h);
          z.setElement(iq, z.getElement(iq) + h);
          d.setElement(ip, dip - h);
          d.setElement(iq, diq + h);
          a.setElement(ip, iq, 0.);

          for (j = 1; j < ip; ++j)
            rotate(a, j - 1, ip, j - 1, iq);
          for (j = ip + 2; j < iq; ++j)
            rotate(a, ip, j - 1, iq, j - 1);
          for (j = iq + 1; j < _width; ++j)
            rotate(a, ip, j, j, iq);
          for (j = 0; j < _width; ++j)
            rotate(v, j, ip, j, iq);
          ++nrot;
        }
      }
    for (ip = 0; ip < _width; ++ip)
    {
      b.setElement(ip, b.getElement(ip) + z.getElement(ip));
      d.setElement(ip, b.getElement(ip));
      z.setElement(ip, 0.);
    }
    }
#ifdef _DEBUG
  _wassert(_CRT_WIDE("Too many iterations in the rotate jacobi"), _CRT_WIDE(__FILE__), (unsigned)(__LINE__));
#endif
  return *this;
}

Matrix& Matrix::rotateJacobi(Matrix& v, Matrix& d) noexcept
{
  *this = rotatedJacobi(v, d);

  return *this;
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
