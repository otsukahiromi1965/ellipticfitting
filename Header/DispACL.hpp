#pragma once

#include <boost/math/tools/roots.hpp> // bisect
extern int LCN;
void disp_ACL_(
    const REAL *para)
{
  if (!rooty)
    return;
  using namespace std;
  using namespace boost::math::tools;
  cout << "begin{disp_ACL}\n";
  const REAL_mp A = para[0];
  const REAL_mp k = para[1];
  const REAL_mp b = para[2];
  const REAL_mp k1 = sqrt(myONE - k * k);
  const REAL_mp K0 = boost::math::ellint_1(k);
  const REAL_mp K1 = boost::math::ellint_1(k1);
  cout << " A = " << A
       << " k = " << k
       << " b = " << b
       << " K0= " << K0
       << " K1= " << K1
       << "\n";
  const REAL_mp phi[4] = {(2 - LCN - b) / LCN, (2 - LCN + b) / LCN, -b / LCN, b / LCN};
  COMPLEX_FUNC_of_REAL_mp v_[4];
  for (int i = 0; i < 4; i++)
    v_[i] = [=](REAL_mp v)
    {
      return K0 + myI * (v + K1 * phi[i]);
    };
  auto sn = complex_sn(k);
  COMPLEX_FUNC_of_REAL_mp complex_f;
  REAL_FUNC_of_REAL_mp f;
  REAL_mp theta;
  REAL_mp s;
  REAL_mp c;
  pair<REAL_mp, REAL_mp> sp;
  REAL_mp vs;
  COMPLEX_mp acl;
  eps_tolerance<REAL_mp> tol(numeric_limits<REAL_mp>::digits);
  vs = 0; // v value for theta == 0
  for (int i = 0; i < 360 / LCN; i++)
  {
    theta = i * 2 * myPI / 360;
    s = sin(theta);
    c = cos(theta - (LCN == 4 ? 0 : myPI / LCN));
    complex_f = [=](REAL_mp v)
    {
      return s * sn(v_[0](v) + v_[1](v)) * sinh(log(k * sn(v_[0](v)) * sn(v_[1](v)))) +
             c * sn(v_[2](v) + v_[3](v)) * sinh(log(k * sn(v_[2](v)) * sn(v_[3](v))));
    };
    f = [=](REAL_mp v)
    {
      return complex_f(v).imag();
    };
    sp = bisect(f, vs - K1 / 9, vs + K1 / 9, tol);
    vs = (sp.first + sp.second) / 2;
    acl = s * log(k * sn(v_[0](vs)) * sn(v_[1](vs))) +
          c * log(k * sn(v_[2](vs)) * sn(v_[3](vs)));
    cout << setw(18) << theta
         << setw(18) << vs
         << setw(18) << abs(complex_f(vs))
         << setw(18) << -sin(2 * myPI / LCN) / acl.real()
         << "\n";
  }
  cout << "end{disp_ACL}\n";
}
