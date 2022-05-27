#include "math.h"

#include <iostream>
#include <ctime> 

int main()
{
  srand(time(nullptr));
  clock_t start = clock();
  const size_t N = 10;
  math::MatrixCSR m = math::MatrixCSR::hilbert(N);
  /*math::MatrixCSR m(N);
  m.reserve(N * N);
  for (size_t i = 0; i < N * N; ++i)
    m.setElement(i, (double)rand() / RAND_MAX * 20 - 10);*/
  math::MatrixCSR mmm1, mmm2;
  std::cout << m.condLU() << std::endl;
  math::MatrixCSR mmm3 = m.rotatedJacobi(mmm1, mmm2);
  std::cout << mmm1 << std::endl;
  std::cout << mmm2 << std::endl;
  std::cout << mmm3 << std::endl;
  /*math::MatrixCSR mInv = m.inverseLU();
  math::MatrixCSR iden = m * mInv;
  std::cout << iden;*/
  clock_t end = clock();
  std::cout << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

  return 0;
}
