#pragma once

// Select fitting samples in (C_min, C_max) and display samples above C_min / 10.

void make_annule_lattice(
    char **argv,
    TWO_DIM *annule,
    TWO_DIM *lattice)
{
  using namespace std;
  LCN = parse_LCN(argv[1]);
  REAL C_MAX = parse_real(argv[3]);
  REAL C_MIN = parse_real(argv[4]);
  REAL C_MIN10 = C_MIN / 10; // cutoff for lattice
  if (!(C_MIN < C_MAX && 0 < C_MIN && C_MAX < 1))
    throw invalid_argument("Expected 0 < C_min < C_max < 1");

  int NX;
  const int NX_MAX = 1000;       // for static array; should be larger than the maximum i, j in datafile
  static REAL c[NX_MAX][NX_MAX]; // Keep the large array off the stack.
  static REAL e[NX_MAX][NX_MAX]; // Keep the large array off the stack.

  if (rooty)
  {
    int i0, j0;
    REAL c0, e0;
    ifstream ifs(argv[2], ios::in);
    if (!ifs)
      throw runtime_error(string("Cannot open input: ") + argv[2]);
    string str;
    NX = 0;
    int line = 0;
    int rows = 0;
    while (getline(ifs, str))
    {
      ++line;
      if (str.find_first_not_of(" \t\r") == string::npos)
        continue;
      stringstream ss(str);
      string extra;
      if (!(ss >> i0 >> j0 >> c0 >> e0) || (ss >> extra) ||
          i0 < 0 || j0 < 0 || i0 >= NX_MAX || j0 >= NX_MAX ||
          !isfinite(c0) || !isfinite(e0) || e0 < 0) // eo can be zero but not negative.
        throw runtime_error("Invalid input row " + to_string(line) +
                            ": expected i j correlation positive_error; indices 0..999");
      ++rows;
      c[i0][j0] = c0;
      e[i0][j0] = e0;
      NX = std::max(NX, std::max(i0, j0));
    }
    if (ifs.bad())
      throw runtime_error("Failed to read input");
    if (rows == 0)
      throw runtime_error("Input is empty");
    NX++;
    ifs.close();
    if (LCN == 4 || LCN == 6)
      for (int i = 0; i < NX; i++) // Reflection across the diagonal.
        for (int j = i + 1; j < NX; j++)
        {
          c0 = c[i][j] + c[j][i];
          e0 = e[i][j] + e[j][i];
          c[i][j] = c0 / 2;
          e[i][j] = e0 / 2;
          c[j][i] = c0 / 2;
          e[j][i] = e0 / 2;
        }
    if (LCN == 6) // Additional symmetry for the triangular lattice.
      for (int i = 0; i < NX; i++)
        for (int j = i + 1; j < NX; j++)
        {
          c0 = c[i][j] + c[j][j - i];
          e0 = e[i][j] + e[j][j - i];
          c[i][j] = c0 / 2;
          e[i][j] = e0 / 2;
          c[j][j - i] = c0 / 2;
          e[j][j - i] = e0 / 2;
        }
  }
  MPI_Bcast(&NX, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(c, NX_MAX * NX_MAX, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(e, NX_MAX * NX_MAX, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);

  lattice->push_back({0, 0, c[0][0], e[0][0]});
  for (int i = 0; i < NX; i++)
    for (int j = 1; j < NX; j++)
    {
      if (C_MIN10 < c[i][j])
        lattice->push_back({i, j, c[i][j], e[i][j]});
      if (C_MIN < c[i][j] && c[i][j] < C_MAX)
        annule->push_back({i, j, c[i][j], e[i][j]});
    }

  if (annule->size() < 3)
    throw runtime_error("At least three fitting samples are required; check cutoffs");
  for (const auto &sample : *annule)
    if (sample.e <= 0)
      throw runtime_error("Fitting samples must have positive errors");

  if (rooty)
    cout << " LCN = " << string(argv[1])
         << "\n DataFile: " << string(argv[2])
         << "\n C_max = " << string(argv[3])
         << "\n C_min = " << string(argv[4])
         << "\n annule.size = " << annule->size()
         << "\n lattice.size = " << lattice->size()
         << endl;
}
