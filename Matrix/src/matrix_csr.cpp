#include "matrix_csr.h"
#include "matrix.h"

#include <cstring>
#include <cassert>

namespace math
{

MatrixCSR::MatrixCSR() noexcept : _width(0), _height(0), _elementCount(0),
  _reservedSize(0), _reservedHeight(0), _data(nullptr), _jptr(nullptr), _iptr(nullptr)
{

}

MatrixCSR::MatrixCSR(const MatrixCSR& matr) noexcept : _width(matr._width), _height(matr._height),
  _elementCount(matr._elementCount), _reservedSize(matr._elementCount), _reservedHeight(matr._height),
  _data(matr._elementCount ? new double[matr._elementCount] : nullptr),
  _jptr(matr._elementCount ? new size_t[matr._elementCount] : nullptr),
  _iptr(matr._height ? new size_t[matr._height + 1] : nullptr)
{
  memcpy(_data, matr._data, sizeof(double) * _elementCount);
  memcpy(_jptr, matr._jptr, sizeof(size_t) * _elementCount);
  if (_iptr)
    memcpy(_iptr, matr._iptr, sizeof(size_t) * (_height + 1));
}

MatrixCSR::MatrixCSR(MatrixCSR&& matr) noexcept : _width(matr._width), _height(matr._height),
  _elementCount(matr._elementCount), _reservedSize(matr._reservedSize), _reservedHeight(matr._reservedHeight),
  _data(matr._data), _jptr(matr._jptr), _iptr(matr._iptr)
{
  matr._width = 0;
  matr._height = 0;
  matr._elementCount = 0;
  matr._reservedSize = 0;
  matr._reservedHeight = 0;

  matr._data = nullptr;
  matr._jptr = nullptr;
  matr._iptr = nullptr;
}

MatrixCSR::MatrixCSR(const Matrix& matr) noexcept : _width(matr.getWidth()), _height(matr.getHeight()),
  _elementCount(0), _reservedSize(0), _reservedHeight(matr.getHeight()), _data(nullptr), _jptr(nullptr), _iptr(nullptr)
{
  size_t matrSize = _width * _height;
  size_t countNotZero = 0;
  for (size_t num = 0; num < matrSize; ++num)
    countNotZero += !!matr[num];
  _elementCount = countNotZero;
  _reservedSize = countNotZero;

  if (!countNotZero)
    return;

  _data = new double[countNotZero];
  _jptr = new size_t[countNotZero];
  _iptr = new size_t[_height + 1];
  _iptr[0] = 0;

  for (size_t y = 1; y < _height + 1; ++y)
    _iptr[y] = countNotZero;

  size_t num = 0;
  for (size_t y = 0; y < _height && num < countNotZero; ++y)
  {
    size_t lineCount = 0;
    for (size_t x = 0; x < matr.getWidth() && num < countNotZero; ++x)
    {
      double elem = matr.getElement(y, x);
      if (elem)
      {
        _data[num] = elem;
        _jptr[num] = x;
        num++;
        ++lineCount;
      }
    }
    _iptr[y + 1] = _iptr[y] + lineCount;
  }
}

MatrixCSR::MatrixCSR(size_t width, size_t height) noexcept : _width(width), _height(height),
  _elementCount(0), _reservedSize(0), _reservedHeight(height), _data(nullptr), _jptr(nullptr),
  _iptr(height ? new size_t[height + 1] : nullptr)
{
  if (_iptr)
    memset(_iptr, 0, sizeof(size_t) * (height + 1));
}

MatrixCSR::MatrixCSR(size_t size) noexcept : _width(size), _height(size), _elementCount(0),
  _reservedSize(0), _reservedHeight(size), _data(nullptr), _jptr(nullptr), _iptr(size ? new size_t[size + 1] : nullptr)
{
  if (_iptr)
    memset(_iptr, 0, sizeof(size_t) * (size + 1));
}

MatrixCSR::~MatrixCSR() noexcept
{
  delete[] _data;
  delete[] _jptr;
  delete[] _iptr;
}

MatrixCSR& MatrixCSR::operator=(const MatrixCSR& matr) noexcept
{
  if (this == &matr)
    return *this;

  if (matr._elementCount)
  {
    if (matr._elementCount > _reservedSize)
    {
      double* newData = new double[matr._elementCount];
      size_t* newJptr = new size_t[matr._elementCount];

      memcpy(newData, matr._data, sizeof(double) * matr._elementCount);
      memcpy(newJptr, matr._jptr, sizeof(size_t) * matr._elementCount);

      delete[] _data;
      delete[] _jptr;

      _data = newData;
      _jptr = newJptr;
      _reservedSize = matr._elementCount;
    }
    else
    {
      memcpy(_data, matr._data, sizeof(double) * matr._elementCount);
      memcpy(_jptr, matr._jptr, sizeof(size_t) * matr._elementCount);
    }
  }

  if (_reservedHeight < matr._height)
  {
    size_t* newIptr = new size_t[matr._height + 1];
    memcpy(newIptr, matr._iptr, sizeof(size_t) * (matr._height + 1));
    delete[] _iptr;
    _iptr = newIptr;
    _reservedHeight = matr._height;
  }
  else
  {
    memcpy(_iptr, matr._iptr, sizeof(size_t) * (matr._height + 1));
  }

  _width = matr._width;
  _height = matr._height;
  _elementCount = matr._elementCount;

  return *this;
}

MatrixCSR& MatrixCSR::operator=(MatrixCSR&& matr) noexcept
{
  if (this == &matr)
    return *this;

  _width = matr._width;
  _height = matr._height;
  _elementCount = matr._elementCount;
  _reservedSize = matr._reservedSize;
  _reservedHeight = matr._reservedHeight;

  matr._width = 0;
  matr._height = 0;
  matr._elementCount = 0;
  matr._reservedSize = 0;
  matr._reservedHeight = 0;

  delete[] _data;
  delete[] _jptr;
  delete[] _iptr;

  _data = matr._data;
  _jptr = matr._jptr;
  _iptr = matr._iptr;

  matr._data = nullptr;
  matr._jptr = nullptr;
  matr._iptr = nullptr;

  return *this;
}

MatrixCSR& MatrixCSR::operator=(const Matrix& matr) noexcept
{
  size_t matrSize = matr.getWidth() * matr.getHeight();
  size_t countNotZero = 0;
  for (size_t num = 0; num < matrSize; ++num)
    countNotZero += !!matr[num];

  if (matr.getHeight() > _reservedHeight)
  {
    size_t* newIptr = new size_t[matr.getHeight() + 1];
    delete[] _iptr;
    _iptr = newIptr;

    _reservedHeight = matr.getHeight();
  }

  if (countNotZero != _elementCount)
  {
    if (countNotZero)
    {
      if (countNotZero > _reservedSize)
      {
        _reservedSize = countNotZero;

        delete[] _data;
        delete[] _jptr;

        _data = new double[countNotZero];
        _jptr = new size_t[countNotZero];
      }
    }
    else
    {
      delete[] _data;
      delete[] _jptr;

      _data = nullptr;
      _jptr = nullptr;
    }
    _elementCount = countNotZero;
  }

  if (countNotZero)
  {
    _iptr[0] = 0;

    for (size_t y = 1; y < matr.getHeight() + 1; ++y)
      _iptr[y] = countNotZero;

    size_t num = 0;
    for (size_t y = 0; y < matr.getHeight() && num < countNotZero; ++y)
    {
      size_t lineCount = 0;
      for (size_t x = 0; x < matr.getWidth() && num < countNotZero; ++x)
      {
        double elem = matr.getElement(y, x);
        if (elem)
        {
          _data[num] = elem;
          _jptr[num] = x;
          ++num;
          ++lineCount;
        }
      }
      _iptr[y + 1] = _iptr[y] + lineCount;
    }
  }

  _width = matr.getWidth();
  _height = matr.getHeight();

  return *this;
}

MatrixCSR MatrixCSR::operator*(const MatrixCSR& matr) const noexcept
{
  assert(_width == matr._height);
  size_t i, j, k;
  double val;
  MatrixCSR res(matr._width, _height);
  res.reserve(res._width * res._height);
  for (i = 0; i < res._height; ++i)
    if (_iptr[i + 1] - _iptr[i])
      for (j = 0; j < res._width; ++j)
      {
        val = 0;
        for (k = _iptr[i]; k < _iptr[i + 1]; ++k)
          val += _data[k] * matr[_jptr[k] * matr._width + j];
        if (val)
          res.setElement(i, j, val);
      }
  res.shrink_to_fit();
  return res;
}

MatrixCSR& MatrixCSR::operator*=(const MatrixCSR& matr) noexcept
{
  *this = *this * matr;
  return *this;
}

MatrixCSR::operator Matrix() const noexcept
{
  Matrix matr(_width, _height);
  for (size_t num = 0; num < _height * _width; ++num)
    matr[num] = getElement(num);
  return matr;
}

MatrixCSR::operator const Matrix() const noexcept
{
  Matrix matr(_width, _height);
  for (size_t num = 0; num < _height * _width; ++num)
    matr[num] = getElement(num);
  return matr;
}

bool MatrixCSR::binarySearch(size_t l, size_t r, size_t key, size_t& ind) const noexcept
{
  if (r - l == _height)
  {
    ind = l + key;
    return true;
  }
  ++l;
  while (l <= r)
  {
    ind = (l + r) / 2;

    if (_jptr[ind - 1] == key)
    {
      --ind;
      return true;
    }
    if (_jptr[ind - 1] > key)
      r = ind - 1;
    else
      l = ind + 1;
  }
  ind = r;
  return false;
}

double MatrixCSR::getElement(size_t y, size_t x) const noexcept
{
  assert(y < _height);
  assert(x < _width);
  size_t ind;
  if (binarySearch(_iptr[y], _iptr[y + 1], x, ind))
    return _data[ind];
  else
    return 0;
}

double MatrixCSR::getElement(size_t num) const noexcept
{
  assert(num < _width * _height);
  return getElement(num / _width, num % _width);
}

void MatrixCSR::setElement(size_t y, size_t x, double value) noexcept
{
  assert(y < _height);
  assert(x < _width);

  size_t ind;
  bool isFound = binarySearch(_iptr[y], _iptr[y + 1], x, ind);
  if (value)
  {
    if (isFound)
    {
      _data[ind] = value;
      return;
    }
    if (_elementCount >= _reservedSize)
    {
      const size_t newReservedSize = _reservedSize > 1 ? _reservedSize + (_reservedSize >> 1) : 2;
      double* newData = new double[newReservedSize];
      size_t* newJptr = new size_t[newReservedSize];

      memcpy(newData, _data, sizeof(double) * ind);
      memcpy(newJptr, _jptr, sizeof(size_t) * ind);

      newData[ind] = value;
      newJptr[ind] = x;

      memcpy(newData + (ind + 1), _data + ind, sizeof(double) * (_elementCount - ind));
      memcpy(newJptr + (ind + 1), _jptr + ind, sizeof(size_t) * (_elementCount - ind));

      delete[] _data;
      delete[] _jptr;

      _data = newData;
      _jptr = newJptr;
      ++_elementCount;
      _reservedSize = newReservedSize;
    }
    else
    {
      ++_elementCount;
      for (size_t index = _elementCount - 1; index > ind; --index)
      {
        _data[index] = _data[index - 1];
        _jptr[index] = _jptr[index - 1];
      }
      _data[ind] = value;
      _jptr[ind] = x;
    }
    for (size_t index = y + 1; index < _height + 1; ++index)
      ++_iptr[index];
  }
  else
  {
    if (isFound)
    {
      --_elementCount;
      for (size_t index = ind; index < _elementCount; ++index)
      {
        _data[index] = _data[index + 1];
        _jptr[index] = _jptr[index + 1];
      }
      for (size_t index = y + 1; index < _height + 1; ++index)
        --_iptr[index];
      return;
    }
  }
}

void MatrixCSR::setElement(size_t num, double value) noexcept
{
  assert(num < _width * _height);
  setElement(num / _width, num % _width, value);
}

double MatrixCSR::operator[](size_t num) const noexcept
{
  return getElement(num);
}

void MatrixCSR::reserve(size_t capacity) noexcept
{
  if (capacity <= _reservedSize)
    return;

  double* newData = new double[capacity];
  size_t* newJptr = new size_t[capacity];

  memcpy(newData, _data, sizeof(double) * _reservedSize);
  memcpy(newJptr, _jptr, sizeof(size_t) * _reservedSize);

  delete[] _data;
  delete[] _jptr;

  _data = newData;
  _jptr = newJptr;

  _reservedSize = capacity;
}

void MatrixCSR::shrink_to_fit() noexcept
{
  if (_elementCount == _reservedSize)
    return;

  if (_elementCount)
  {
    double* newData = new double[_elementCount];
    size_t* newJptr = new size_t[_elementCount];

    memcpy(newData, _data, sizeof(double) * _elementCount);
    memcpy(newJptr, _jptr, sizeof(size_t) * _elementCount);

    delete[] _data;
    delete[] _jptr;

    _data = newData;
    _jptr = newJptr;
  }
  else
  {
    delete[] _data;
    delete[] _jptr;

    _data = nullptr;
    _jptr = nullptr;
  }

  _reservedSize = _elementCount;
}

size_t MatrixCSR::size() const noexcept
{
  return _elementCount;
}

size_t MatrixCSR::capacity() const noexcept
{
  return _reservedSize;
}

size_t MatrixCSR::getWidth() const noexcept
{
  return _width;
}

size_t MatrixCSR::getHeight() const noexcept
{
  return _height;
}

void MatrixCSR::lu(MatrixCSR& l, MatrixCSR& u) const noexcept
{
  assert(_width == _height);
  l.reserve(_width * _height);
  u.reserve(_width * _height);
  u = *this;
  l = MatrixCSR(_width);

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

MatrixCSR MatrixCSR::solveSystem(const MatrixCSR& system, const MatrixCSR& answers) noexcept
{
  assert(system._width == system._height);
  assert(system._height == answers._height);
  assert(answers._width == 1);

  MatrixCSR l, u;
  system.lu(l, u);
  return solveSystem(l, u, answers);
}

MatrixCSR MatrixCSR::solveSystem(const MatrixCSR& l, const MatrixCSR& u, const MatrixCSR& answers) noexcept
{
  assert(l._width == l._height);
  assert(l._width == u._width);
  assert(u._width == u._height);
  assert(l._height == answers._height);
  assert(answers._width == 1);

  MatrixCSR tmpRes(1, l._width);
  tmpRes.reserve(l._width);
  size_t i, j, k;
  for (i = 0; i < tmpRes.getHeight(); ++i)
  {
    double sum = 0;
    for (j = 0, k = l._iptr[i]; j < i && k < l._iptr[i + 1]; ++j, ++k)
      sum += l._data[k] * tmpRes[j];
    tmpRes.setElement(i, answers[i] - sum);
  }

  MatrixCSR res(1, l._width);
  res.reserve(l._width);
  for (i = res.getHeight(); i > 0; --i)
  {
    double sum = 0;
    for (j = res.getHeight(), k = u._iptr[i]; j > i && k > u._iptr[i - 1]; --j, --k)
      sum += u._data[k - 1] * res[j - 1];
    res.setElement(i - 1, (tmpRes[i - 1] - sum) / u.getElement(i - 1, i - 1));
  }
  return res;
}

MatrixCSR MatrixCSR::inverseLU() const noexcept
{
  assert(_width == _height);
  MatrixCSR inversed(_width, _height);
  MatrixCSR ans(1, _width);
  MatrixCSR l, u;
  lu(l, u);
  size_t i, j;
  inversed.reserve(_width * _height);
  ans.reserve(1);
  ans._data[0] = 1;
  ans._jptr[0] = 0;
  for (i = 0; i < _width; ++i)
  {
    ans._iptr[i] = 0;
    ans._iptr[i + 1] = 1;
    const MatrixCSR solved = solveSystem(l, u, ans);
    for (j = 0; j < _height; ++j)
      inversed.setElement(j, i, solved[j]);
  }
  return inversed;
}

double MatrixCSR::determinantLU() const noexcept
{
  MatrixCSR l, u;
  lu(l, u);
  return determinantLU(u);
}

double MatrixCSR::determinantLU(const MatrixCSR& u) const noexcept
{
  assert(u._width == u._height);
  double res = 1;
  for (size_t i = 0; i < u._width; ++i)
    res *= getElement(i, i);
  return res;
}

double MatrixCSR::condLU() const noexcept
{
  double norm = 0;
  double sum;
  for (size_t i = 0; i < _height; ++i)
  {
    sum = 0;
    for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j)
      sum += abs(_data[j]);
    norm = sum > norm ? sum : norm;
  }
  MatrixCSR inversed = inverseLU();
  double normInv = 0;
  for (size_t i = 0; i < _height; ++i)
  {
    sum = 0;
    for (size_t j = _iptr[i]; j < _iptr[i + 1]; ++j)
      sum += abs(inversed._data[j]);
    normInv = sum > normInv ? sum : normInv;
  }
  return norm * normInv;
}

MatrixCSR& MatrixCSR::transpose() noexcept
{
  *this = transposing();
  return *this;
}

MatrixCSR MatrixCSR::transposing() const noexcept
{
  MatrixCSR result(_height, _width);
  result.reserve(_reservedSize);
  size_t i, j;
  for (i = 0; i < result._height; ++i)
    for (j = 0; j < result._width; ++j)
      result.setElement(i, j, getElement(j, i));
  return result;
}

MatrixCSR MatrixCSR::identity(size_t size) noexcept
{
  MatrixCSR matr(size);
  matr._data = new double[size];
  matr._jptr = new size_t[size];
  for (size_t i = 0; i < size; ++i)
  {
    matr._data[i] = 1;
    matr._jptr[i] = i;
    matr._iptr[i] = i;
  }
  matr._iptr[size] = size;
  matr._elementCount = size;
  matr._reservedSize = size;
  return matr;
}

MatrixCSR MatrixCSR::hilbert(size_t size) noexcept
{
  MatrixCSR matr(size);
  matr._data = new double[size * size];
  matr._jptr = new size_t[size * size];
  matr._iptr[size] = size * size;
  matr._elementCount = size * size;
  matr._reservedSize = size * size;

  size_t i, j;
  for (i = 0; i < size; ++i)
  {
    matr._iptr[i] = i * size;
    for (j = 0; j < size; ++j)
    {
      matr._data[i * size + j] = 1 / static_cast<double>(i + j + 1);
      matr._jptr[i * size + j] = j;
    }
  }
  return matr;
}

std::ostream& operator<<(std::ostream& out, const MatrixCSR& matr) noexcept
{
  out << std::endl;
  for (size_t y = 0; y < matr.getHeight(); ++y)
  {
    for (size_t x = 0; x < matr.getWidth(); ++x)
      out << matr.getElement(y, x) << " ";
    out << std::endl;
  }
  return out;
}

}
