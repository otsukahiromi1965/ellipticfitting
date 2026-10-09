#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <complex>
#include <functional>

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>

#include "Header/MyTypedef.hpp"
#include "Header/InputValidation.hpp"

const REAL_mp myZERO = REAL_mp(0);
const REAL_mp myONE = REAL_mp(1);
const REAL_mp myPI = boost::math::constants::pi<REAL_mp>();
const COMPLEX_mp myI(0, 1);
const int NPARA = 3; // A, k, b
const bool rooty = true;

int LCN; // Lattice coordination number: square = 4, triangular = 6.

#include "Header/SpecialFunction.hpp"
#include "Header/DispACL.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  if (argc == 2 && string(argv[1]) == "--help")
  {
    cout << "Usage: main_acl LCN k b\n";
    return 0;
  }
  try
  {
    if (argc != 4)
      throw invalid_argument("Usage: main_acl LCN k b");
    LCN = parse_LCN(argv[1]);
    // ACL does not depend on the amplitude A.
    const REAL para[NPARA] = {1.0, parse_real(argv[2]), parse_real(argv[3])};
    validate_k(para[1]);
    cout << setprecision(10);
    disp_ACL_(para);
  }
  catch (const exception &error)
  {
    cerr << error.what() << endl;
    return 1;
  }
  return 0;
}
