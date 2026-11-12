#pragma once

#include <complex>
#include <functional>
#include <vector>

#include <boost/multiprecision/cpp_complex.hpp>

using REAL = long double;
using COMPLEX = std::complex<REAL>;
using COMPLEX_mp = boost::multiprecision::cpp_complex_oct;
using REAL_mp = boost::multiprecision::cpp_complex_oct::value_type;
using COMPLEX_FUNC_of_COMPLEX_mp = std::function<COMPLEX_mp(COMPLEX_mp)>;
using COMPLEX_FUNC_of_REAL_mp = std::function<COMPLEX_mp(REAL_mp)>;
using REAL_FUNC_of_REAL_mp = std::function<REAL_mp(REAL_mp)>;

// One lattice sample: coordinates, correlation, and measurement error.
struct MCDATA
{
  int i;
  int j;
  REAL c;
  REAL e;
};

using TWO_DIM = std::vector<MCDATA>;
