#include "math.h"

#include <iostream>
#include <ctime> 

int main()
{
  clock_t start = clock();
  const size_t N = 100;
  math::MatrixCSR m = math::MatrixCSR::hilbert(N);
  /*math::MatrixCSR mInv = m.inverseLU();
  math::MatrixCSR iden = m * mInv;
  std::cout << iden;*/
  std::cout << m.condLU() << std::endl;
  clock_t end = clock();
  std::cout << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

  return 0;
}
