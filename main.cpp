#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>
#include "metod.hpp"
#include "funzioni.cpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>


int main(int argc, char **argv)
{
  TApplication app("app", &argc, argv);
  std::cout << "Avvio myMacro()..." << std::endl;

  gf::prima_prova();
  gf::rigenerazione(100, 1000000, 100);
  gf::bin_smearing(100, 100);
  gf::propagazione_parametri(100, 100);
  std::cout << "myMacro() terminata." << std::endl;
  app.Run();
  return 0;
}
