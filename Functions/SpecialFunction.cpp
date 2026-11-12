#pragma once

COMPLEX_FUNC_of_COMPLEX_mp complex_sn(const REAL_mp k)
{
  const REAL_mp k1 = sqrt(myONE - k * k);
  return [k, k1](COMPLEX_mp z)
  {
    static REAL_mp s0, c0, d0, s1, c1, d1, den;
    s0 = boost::math::jacobi_elliptic(k, z.real(), &c0, &d0);
    s1 = boost::math::jacobi_elliptic(k1, z.imag(), &c1, &d1);
    den = c1 * c1 + k * k * s0 * s0 * s1 * s1;
    const COMPLEX_mp complex_sn(s0 * d1, c0 * d0 * s1 * c1);
    return complex_sn / den;
  };
}

COMPLEX_FUNC_of_REAL_mp complex_sn(const REAL_mp k, const REAL_mp s1, const REAL_mp c1, const REAL_mp d1)
{
  return [=](REAL_mp u)
  {
    static REAL_mp s0, c0, d0, den;
    s0 = boost::math::jacobi_elliptic(k, u, &c0, &d0);
    den = c1 * c1 + k * k * s0 * s0 * s1 * s1;
    const COMPLEX_mp complex_sn(s0 * d1, c0 * d0 * s1 * c1);
    return complex_sn / den;
  };
}

COMPLEX_FUNC_of_REAL_mp complex_cc(const REAL_mp k, const REAL_mp s1, const REAL_mp c1, const REAL_mp d1)
{
  return [=](REAL_mp u)
  {
    static REAL_mp s0, c0, d0, den;
    s0 = boost::math::jacobi_elliptic(k, u, &c0, &d0);
    den = c1 * c1 + k * k * s0 * s0 * s1 * s1;
    const COMPLEX_mp complex_cn(c0 * c1, -s0 * d0 * s1 * d1);
    return (complex_cn / den) * (complex_cn / den);
  };
}

COMPLEX_FUNC_of_REAL_mp complex_cd(const REAL_mp k, const REAL_mp s1, const REAL_mp c1, const REAL_mp d1)
{
  return [=](REAL_mp u)
  {
    static REAL_mp s0, c0, d0, den;
    s0 = boost::math::jacobi_elliptic(k, u, &c0, &d0);
    den = c1 * c1 + k * k * s0 * s0 * s1 * s1;
    const COMPLEX_mp complex_sn(s0 * d1, c0 * d0 * s1 * c1);
    const COMPLEX_mp complex_cn(c0 * c1, -s0 * d0 * s1 * d1);
    const COMPLEX_mp complex_dn(d0 * c1 * d1, -k * k * s0 * c0 * s1);
    return (complex_cn / den) * (complex_dn / den) / (complex_sn / den);
  };
}

COMPLEX_FUNC_of_REAL_mp complex_dz(const REAL_mp k, const REAL_mp phi, const REAL_mp s1, const REAL_mp c1, const REAL_mp d1)
{
  const REAL_mp sc1 = s1 / c1;
  const REAL_mp k1 = sqrt(myONE - k * k);
  const REAL_mp z1 = boost::math::jacobi_zeta(k1, asin(s1)) + phi - d1 * sc1;
  return [=](REAL_mp u)
  {
    static REAL_mp s0, c0, d0, den, z0;
    s0 = boost::math::jacobi_elliptic(k, u, &c0, &d0);
    z0 = boost::math::jacobi_zeta(k, asin(s0)); // Z(k,x) = Z(x|k) = Z(k,phi=asin(sn(k,x)))
    den = c1 * c1 + k * k * s0 * s0 * s1 * s1;
    const COMPLEX_mp complex_sn(s0 * d1, c0 * d0 * s1 * c1);
    const COMPLEX_mp complex_cn(c0 * c1, -s0 * d0 * s1 * d1);
    const COMPLEX_mp complex_dn(d0 * c1 * d1, -k * k * s0 * c0 * s1);
    const COMPLEX_mp complex_zn = z0 - myI * (z1 + k * k * s0 * sc1 * (complex_sn / den));
    return ((complex_cn / den) * (complex_dn / den) / (complex_sn / den)) * complex_zn;
  };
}
