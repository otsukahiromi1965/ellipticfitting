#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <mpi.h>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>

#include "Header/InputValidation.hpp"
#include "Header/MyTypedef.hpp"

const REAL_mp myZERO = REAL_mp(0);
const REAL_mp myONE = REAL_mp(1);
const REAL_mp myPI = boost::math::constants::pi<REAL_mp>();
const COMPLEX_mp myI(0, 1);
const int NPARA = 3; // A, k, b
const int root = 0;

int LCN;     // LatticeCoordinationNumber: square = 4, triangular = 6.
int nu_proc; // Number of MPI processes.
int my_proc; // MPI rank of this process.
bool rooty;  // True if this process is the root process.

#include "Header/SpecialFunction.hpp"
#include "Header/AnnularRegion.hpp"
#include "Header/Gradient.hpp"
#include "Header/DispDeviation.hpp"
#include "Header/DispACL.hpp"

void calculate_residuals(
    const TWO_DIM &system,
    const REAL *fit,
    REAL *res)
{
  if (!rooty)
    return;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
    res[samp] = (system[samp].c - fit[samp]) / system[samp].e;
}

void calculate_mean_squared_residual(
    const TWO_DIM &system,
    const REAL *res,
    REAL *chi)
{
  if (!rooty)
    return;
  REAL x = 0;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
    x += pow(res[samp], 2);
  *chi = x / static_cast<int>(system.size());
}

void calculate_parameter_step(
    const TWO_DIM &system,
    const REAL *res,
    const std::array<std::vector<REAL>, NPARA> &gradient,
    const REAL lambda,
    REAL *diff)
{
  if (!rooty)
    return;
  using namespace Eigen;
  using MatrixREAL = Matrix<REAL, NPARA, NPARA, RowMajor>;
  using VectorREAL = Matrix<REAL, NPARA, 1>;
  MatrixREAL A = MatrixREAL::Zero(NPARA, NPARA);
  VectorREAL v = VectorREAL::Zero(NPARA);
  VectorREAL g;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
  {
    for (int i = 0; i < NPARA; i++)
      g(i) = gradient[i][samp] / system[samp].e;
    A += g * g.transpose();
    v += g * res[samp];
  }
  A += lambda * A.diagonal().asDiagonal();
  Map<VectorREAL>(diff, NPARA) = A.fullPivLu().solve(v);
}

void gather_sample_values(
    const TWO_DIM &system,
    const REAL *local_values,
    REAL *gathered_values)
{
  // Each sample belongs to exactly one rank. Reduce avoids blocking self-sends.
  std::vector<REAL> local(system.size(), 0);
  for (size_t samp = my_proc; samp < system.size(); samp += nu_proc)
    local[samp] = local_values[samp];
  MPI_Reduce(local.data(), gathered_values, static_cast<int>(system.size()),
             MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  rooty = (my_proc == root);

  using namespace std;
  if (argc == 2 && string(argv[1]) == "--help")
  {
    if (rooty)
      cout << "Usage: main LCN datafile C_max C_min A k b\n"
           << "LCN: 4 or 6; 0 < C_min < C_max < 1; 0 < k < 1\n";
    MPI_Finalize();
    return 0;
  }
  try
  {
    if (argc != 8)
      throw invalid_argument("Usage: main LCN datafile C_max C_min A k b");
    LCN = parse_LCN(argv[1]);
    const REAL cmax = parse_real(argv[3]);
    const REAL cmin = parse_real(argv[4]);
    if (!(0 < cmin && cmin < cmax && cmax < 1))
      throw invalid_argument("Expected 0 < C_min < C_max < 1");
    const int max_iterations = 100;
    const REAL relative_tolerance = 1e-6;
    REAL chi[2] = {}; // Residual scores for the two parameter slots.
    REAL diff[NPARA];
    REAL lambda = 100;
    const REAL initial_parameters[NPARA] = {
        parse_real(argv[5]), parse_real(argv[6]), parse_real(argv[7])};
    REAL para[2][NPARA]; // Two slots for candidate and accepted parameter sets.
    for (auto &parameters : para)
      std::copy_n(initial_parameters, NPARA, parameters);
    validate_k(initial_parameters[1]);
    // The candidate and accepted parameter sets alternate between slots 0 and 1.
    bool accepted = false;
    TWO_DIM annule;
    TWO_DIM lattice;
    make_annule_lattice(argv, &annule, &lattice);

    std::vector<REAL> local_fit(annule.size());
    std::vector<REAL> gathered_fit(annule.size());
    std::vector<REAL> local_lattice_fit(lattice.size());
    std::vector<REAL> gathered_lattice_fit(lattice.size());
    std::array<std::vector<REAL>, NPARA> local_gradient;
    std::array<std::vector<REAL>, NPARA> gathered_gradient;
    std::array<std::vector<REAL>, 2> residuals;
    for (auto &gradient : local_gradient)
      gradient.resize(annule.size());
    for (auto &gradient : gathered_gradient)
      gradient.resize(annule.size());
    for (auto &residual : residuals)
      residual.resize(annule.size());

    if (rooty)
    {
      cout << "begin{loop}\n";
      cout << " Loop; "
           << " A; "
           << " k; "
           << " b; "
           << " X; "
           << " Lambda\n";
    }

    cout << setprecision(10);

    int loop = 0;
    while (loop < max_iterations)
    {
      calculate_f(annule, para[!accepted], local_fit.data());
      gather_sample_values(annule, local_fit.data(), gathered_fit.data());
      calculate_residuals(annule, gathered_fit.data(), residuals[!accepted].data());
      calculate_mean_squared_residual(annule, residuals[!accepted].data(), &chi[!accepted]);
      MPI_Bcast(chi, 2, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);

      if (loop == 0 || chi[!accepted] < chi[accepted])
      {
        accepted = !accepted;

        if (rooty)
        {
          cout.setf(std::ios::left);
          cout << setw(3) << loop << " "
               << setw(14) << para[accepted][0] << " "
               << setw(14) << para[accepted][1] << " "
               << setw(14) << para[accepted][2] << " "
               << setw(15) << chi[accepted] << " "
               << lambda << " "
               << endl;
        }

        if (chi[accepted] == 0 ||
            (loop != 0 && (chi[!accepted] - chi[accepted]) / chi[accepted] < relative_tolerance))
          break;

        calculate_gradient(annule, para[accepted], local_fit.data(),
                           local_gradient[0].data(), local_gradient[1].data(), local_gradient[2].data());
        for (int i = 0; i < NPARA; i++)
          gather_sample_values(annule, local_gradient[i].data(), gathered_gradient[i].data());

        lambda /= 10;
      }
      else
        lambda *= 10;

      calculate_parameter_step(annule, residuals[accepted].data(), gathered_gradient, lambda, diff);
      MPI_Bcast(diff, NPARA, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
      for (int i = 0; i < NPARA; i++)
        para[!accepted][i] = para[accepted][i] + diff[i];
      for (int i = 0; i < NPARA; ++i)
        if (!std::isfinite(para[!accepted][i]))
          throw runtime_error("Non-finite fit parameter; try different initial values");
      validate_k(para[!accepted][1]);

      ++loop;
    }
    if (rooty)
      cout << "end{loop}\n";

    calculate_f(lattice, para[accepted], local_lattice_fit.data());
    gather_sample_values(lattice, local_lattice_fit.data(), gathered_lattice_fit.data());
    disp_deviation_(loop, lattice, gathered_lattice_fit.data());
    disp_ACL_(para[accepted]);
  }
  catch (const std::exception &error)
  {
    cerr << "Rank " << my_proc << ": " << error.what() << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  MPI_Finalize();
  return 0;
}
