#pragma once

void disp_deviation_(
    const int loop, // for display
    const TWO_DIM &system,
    const REAL *fit)
{
  if (!rooty)
    return;
  using namespace std;
  cout << "begin{disp_deviation_" << loop << "}" << "\n";
  static REAL data;
  static REAL form;
  for (int samp = 0; samp < static_cast<int>(system.size()); samp++)
  {
    form = fit[samp];
    data = system[samp].c;
    cout << setw(4) << system[samp].i
         << setw(4) << system[samp].j
         << setw(18) << data
         << setw(18) << form
         << setw(18) << ((form - data) / data)
         << "\n";
  }
  cout << "end{disp_deviation_" << loop << "}\n";
}
