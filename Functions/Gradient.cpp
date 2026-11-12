#pragma once

#include <boost/math/quadrature/tanh_sinh.hpp>
boost::math::quadrature::tanh_sinh<REAL_mp> integrator;
#include <iostream>

// Integer powers avoid unnecessary complex logarithms.
template <typename T>
T my_pow(T base, int exponent)
{
  T result = 1;
  while (exponent > 0)
  {
    if (exponent & 1)
      result *= base;
    base *= base;
    exponent >>= 1;
  }
  return result;
}

// Half-period integration for the model and its gradient.

void calculate_f(
    const TWO_DIM &system,
    const REAL *para,
    REAL *fit)
{
  using namespace std;
  const REAL_mp A = para[0];
  const REAL_mp k = para[1];
  const REAL_mp b = para[2];
  const REAL_mp k1 = sqrt(myONE - k * k);
  const REAL_mp K0 = boost::math::ellint_1(k);
  const REAL_mp K1 = boost::math::ellint_1(k1);
  const REAL_mp amp = A * sqrt(k1) / myPI;
  const REAL_mp phi[4] = {(2 - LCN - b) / LCN, (2 - LCN + b) / LCN, -b / LCN, b / LCN};

  static COMPLEX_FUNC_of_REAL_mp sn_[4];
  for (int i = 0; i < 4; i++)
  {
    REAL_mp s1, c1, d1;
    s1 = boost::math::jacobi_elliptic(k1, K1 * phi[i], &c1, &d1);
    sn_[i] = complex_sn(k, s1, c1, d1);
  }

  static COMPLEX_FUNC_of_REAL_mp xbase;
  static COMPLEX_FUNC_of_REAL_mp ybase;
  xbase = [=](REAL_mp u)
  {
    return (k * sn_[0](u) * sn_[1](u));
  };
  ybase = [=](REAL_mp u)
  {
    return (k * sn_[2](u) * sn_[3](u));
  };

  static REAL_mp v;
  static COMPLEX_mp complex_f;
  static REAL_FUNC_of_REAL_mp f;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
  {
    if (my_proc != samp % nu_proc)
      continue;
    f = [=](REAL_mp u)
    {
      v = K0 - u;
      complex_f =
          my_pow(xbase(u), system[samp].i) *
          my_pow(ybase(u), system[samp].j) +
          my_pow(xbase(v), system[samp].i) *
          my_pow(ybase(v), system[samp].j);
      return complex_f.real();
    };
    fit[samp] = (REAL)(amp * 2 * integrator.integrate(f, myZERO, K0 / 2));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void calculate_gradient(
    const TWO_DIM &system,
    const REAL *para,
    const REAL *fit,
    REAL *DfDA,
    REAL *DfDk,
    REAL *DfDb)
{
  using namespace std;
  const REAL_mp A = para[0];
  const REAL_mp k = para[1];
  const REAL_mp b = para[2];
  const REAL_mp m = k * k;
  const REAL_mp m1 = myONE - m;
  const REAL_mp k1 = sqrt(m1);
  const REAL_mp K0 = boost::math::ellint_1(k);
  const REAL_mp K1 = boost::math::ellint_1(k1);
  const REAL_mp E0 = boost::math::ellint_2(k);
  const REAL_mp amp = A * sqrt(k1) / myPI;
  const REAL_mp phi[4] = {(2 - LCN - b) / LCN, (2 - LCN + b) / LCN, -b / LCN, b / LCN};
  const REAL_mp C0 = myPI / 2 / K0;
  const REAL_mp C1 = myPI / 2 / m1 / K0;
  const REAL_mp C2 = -(E0 - m1 * K0) / m1 / K0;
  const REAL_mp C3 = m / m1;
  const REAL_mp C4 = -myONE / m1;
  const COMPLEX_mp myPL = myI * (-phi[0] - phi[1]) / 2;
  const COMPLEX_mp myMI = myI * (-phi[0] + phi[1]) / 2;

  static COMPLEX_FUNC_of_REAL_mp sn_[4];
  static COMPLEX_FUNC_of_REAL_mp cc_[4];
  static COMPLEX_FUNC_of_REAL_mp cd_[4];
  static COMPLEX_FUNC_of_REAL_mp dz_[4];
  for (int i = 0; i < 4; i++)
  {
    REAL_mp s1, c1, d1;
    s1 = boost::math::jacobi_elliptic(k1, K1 * phi[i], &c1, &d1);
    sn_[i] = complex_sn(k, s1, c1, d1);
    cc_[i] = complex_cc(k, s1, c1, d1);
    cd_[i] = complex_cd(k, s1, c1, d1);
    dz_[i] = complex_dz(k, C0 * phi[i], s1, c1, d1);
  }

  static COMPLEX_FUNC_of_REAL_mp xbase;
  static COMPLEX_FUNC_of_REAL_mp ybase;
  static COMPLEX_FUNC_of_REAL_mp xcoef;
  static COMPLEX_FUNC_of_REAL_mp ycoef;
  xbase = [=](REAL_mp u)
  {
    return (k * sn_[0](u) * sn_[1](u));
  };
  ybase = [=](REAL_mp u)
  {
    return (k * sn_[2](u) * sn_[3](u));
  };
  xcoef = [=](REAL_mp u)
  {
    return (C1 * (cd_[0](u) + cd_[1](u)) * myPL +
            C1 * (cd_[0](u) - cd_[1](u)) * myMI +
            C2 * (cd_[0](u) + cd_[1](u)) * u +
            C3 * (cc_[0](u) + cc_[1](u)) +
            C4 * (dz_[0](u) + dz_[1](u)));
  };
  ycoef = [=](REAL_mp u)
  {
    return (C1 * (cd_[2](u) - cd_[3](u)) * myMI +
            C2 * (cd_[2](u) + cd_[3](u)) * u +
            C3 * (cc_[2](u) + cc_[3](u)) +
            C4 * (dz_[2](u) + dz_[3](u)));
  };

  static REAL_mp x;
  static REAL_mp y;
  static REAL_mp v;
  static REAL_mp DfDk_mp;
  static COMPLEX_mp DfDk_c;
  static COMPLEX_mp complex_f;
  static REAL_FUNC_of_REAL_mp f;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
  {
    if (my_proc != samp % nu_proc)
      continue;
    x = system[samp].i;
    y = system[samp].j;

    DfDA[samp] = fit[samp] / para[0];

    DfDk_c =
        my_pow(xbase(K0), system[samp].i) *
        my_pow(ybase(K0), system[samp].j);
    f = [=](REAL_mp u)
    {
      v = K0 - u;
      complex_f =
          (x * xcoef(u) + y * ycoef(u)) *
              my_pow(xbase(u), system[samp].i) *
              my_pow(ybase(u), system[samp].j) +
          (x * xcoef(v) + y * ycoef(v)) *
              my_pow(xbase(v), system[samp].i) *
              my_pow(ybase(v), system[samp].j);
      return complex_f.real();
    };
    DfDk_mp = -K0 * C2 * DfDk_c.real() + integrator.integrate(f, myZERO, K0 / 2);
    DfDk[samp] = (REAL)((x + y - C3 / 2) / k) * fit[samp] + (REAL)(amp * 2 * DfDk_mp / k);

    if (LCN == 6)
      continue;
    f = [=](REAL_mp u)
    {
      v = K0 - u;
      complex_f =
          (x * (cd_[0](u) - cd_[1](u)) + y * (cd_[2](u) - cd_[3](u))) *
              my_pow(xbase(u), system[samp].i) *
              my_pow(ybase(u), system[samp].j) +
          (x * (cd_[0](v) - cd_[1](v)) + y * (cd_[2](v) - cd_[3](v))) *
              my_pow(xbase(v), system[samp].i) *
              my_pow(ybase(v), system[samp].j);
      return complex_f.imag();
    };
    DfDb[samp] = (REAL)((K1 / LCN) * amp * 2 * integrator.integrate(f, myZERO, K0 / 2));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
