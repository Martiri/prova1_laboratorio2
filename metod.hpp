#ifndef METODO_HPP  
#define METODO_HPP  
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>

namespace gf{
class Function
{
public:
  double k() const;
  double phi() const;
  double b() const;

  TF1 *f;
  TCanvas *c;

  Function();
  

  ~Function();

  TH1D *GH(int N, int nbins);

  double ssqm(TH1D *h, TF1 *f);

  void Draw();

};

void prima_prova();
void rigenerazione(int M, int N, int nbins);
void bin_smearing(int M, int nbins);
void propagazione_parametri(int M, int nbins);
void fit_distribution();
} //namespace gf

#endif
