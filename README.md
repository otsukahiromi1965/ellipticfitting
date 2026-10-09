# EllipticFitting

A C++ program for fitting correlation functions using elliptic functions.
Sample values and gradients are computed in parallel with MPI, and the
Levenberg–Marquardt method estimates the parameters `A`, `k`, and `b`.
Supported lattice coordination numbers are `LCN=4` (square) and `LCN=6`
(triangular). The `main_acl` executable computes the angular dependence of
the correlation length (ACL) from specified parameters without fitting.

## Requirements and build

- A C++ compiler supporting C++23
- CMake 3.16 or later
- Boost headers
- Eigen 3.3 or later
- MPI C++ development libraries and a launcher, such as Open MPI

After installing the dependencies, run these commands from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```

Compiling the multiprecision elliptic functions can take significant time and
memory. Use `--parallel 1` if the build runs out of memory.

## Fitting

```text
mpirun -np <process_count> ./build/main LCN datafile C_max C_min A k b
```

| Argument | Description |
| --- | --- |
| `LCN` | Lattice coordination number: `4` or `6` |
| `datafile` | Path to the input text file |
| `C_max`, `C_min` | Correlation range used for fitting: `0 < C_min < C_max < 1` |
| `A`, `k`, `b` | Initial parameters; all must be finite, with `0 < k < 1` |

Example using the included data:

```sh
mpirun -np 2 ./build/main 6 Data/6/0/NX28_x20 5E-2 1E-2 1.0 0.567 1.0 > result.txt
```

The input file must contain four whitespace-separated columns without a header.
Blank lines are ignored.

```text
i j correlation error
```

The indices `i` and `j` must be integers between `0` and `999`, inclusive.
Correlation values and errors must be finite. Input errors may be zero but
must not be negative; errors for samples selected for fitting must be positive
after symmetrization. If a coordinate appears more than once, the last row is used.

The program assumes coordinate data covering a square region and treats missing
coordinates as zero. Supply data for the target region, including coordinates
referenced during symmetrization. Correlations and errors are averaged according
to the lattice symmetry, then points with `j >= 1` and
`C_min < correlation < C_max` are selected for fitting. At least three fitting
points are required. The implementation uses fixed arrays with 1,000 entries
along each axis; indices outside this range are rejected.

Results are written to standard output:

- `begin{loop}`: iteration number, `A`, `k`, `b`, mean squared normalized
  residual `X`, and damping coefficient `Lambda`.
- `begin{disp_deviation_...}`: coordinates, input correlations, fitted values,
  and relative deviations.
- `begin{disp_ACL}`: parameters, followed by angles in radians, the root-solved
  value of `v`, the absolute function residual, and the correlation length.

The maximum number of iterations is 100, and the relative improvement threshold
is `1e-6`. Depending on the initial parameters, fitting may fail to converge,
or elliptic-function evaluation or ACL root finding may fail. The program exits
with an error if an updated `k` is outside its valid range. Adjust the initial
parameters or fitting range if necessary.

The repository also includes `go0.sh`, which accepts optional input and output
paths:

```sh
./go0.sh path/to/input.txt path/to/result.txt
```

This script uses the host file `MPI_Local_Host`, sets `LCN=6` and fixed fitting
parameters, and launches the MPI job in the background. Edit the script to match
your host configuration and fitting settings before using it.

## Computing ACL only

```sh
./build/main_acl 6 0.567 1.0
```

The arguments are `LCN k b`. The amplitude `A` is not used in ACL calculations.
Both executables display usage information when invoked with `--help`.

## Publishing and licensing

`Data/` contains research data and existing calculation results. Select only
data suitable for public release when publishing the repository.
The `.gitignore` file excludes matching output files, build artifacts, and
VS Code settings, but it does not affect files already tracked by Git.
Check local host files separately before publishing.

No license is currently specified for this repository. The rights holders of
the code and data should determine and document the terms for reuse and
redistribution.
