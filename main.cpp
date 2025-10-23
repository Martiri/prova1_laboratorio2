#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

/*TRandom3* rng = new TRandom3();
rng->SetSeed();
*/
class Function
{
public:
  double k{5.2};
  double phi{1.8};
  double b{0.2};

  TF1 *f;
  TCanvas *c;

  Function()
  {
    f = new TF1("f", "pow(cos([0]*x + [1]),2) + [2]", 0., 4.);
    f->SetParNames("k", "phi", "b");
    f->SetParameters(k, phi, b);
    c = nullptr;
  }

  ~Function()
  {
    delete f;
    if (c)
      delete c;
  }

  TH1D *GH(int N, int nbins)
  {
    // Creiamo un nome unico per ogni istogramma
    static int counter = 0;
    std::string histName = "h_" + std::to_string(counter++);

    TH1D *h = new TH1D(histName.c_str(), "Distribuzione generata", nbins, 0., 4.);
    for (int i = 0; i < N; ++i)
    {
      double x = f->GetRandom();
      h->Fill(x);
    }

    return h;
  }
  double ssqm(TH1D *h, TF1 *f)
  { // mi fa la deviazione standard
    int nbins = h->GetNbinsX();
    std::vector<double> sqm;
    sqm.reserve(nbins);
    for (int i = 1; i <= nbins; ++i)
    {
      double y_h = h->GetBinContent(i);
      double y_f = f->Eval(h->GetBinCenter(i));
      double d = y_h - y_f;
      sqm.push_back(d * d);
    }

    double sum = std::accumulate(sqm.begin(), sqm.end(), 0.);
    double dev = std::sqrt(sum / nbins);
    return dev;
  }

  void Draw()
  {
    c = new TCanvas("c", "Grafico cos^2(kx+phi)+b", 800, 600);
    f->SetTitle("f(x) = cos^{2}(k x + #varphi) + b; x; f(x)");
    f->Draw();
    c->Update();
  }
};

void prima_prova()
{
  static Function
      myfun; // serve per evitare che venga distrutta alla fine della funzione
  myfun.Draw();
  myfun.c->Update();
  TH1D *h = myfun.GH(1000000, 100);
  h->SetLineColor(kBlue);

  // Normalizza l'istogramma all'area della funzione
  double x_min = h->GetXaxis()->GetXmin();
  double x_max = h->GetXaxis()->GetXmax();
  double fIntegral = myfun.f->Integral(x_min, x_max);
  double hIntegral = h->Integral("width");
  if (hIntegral > 0.0 && fIntegral > 0.0)
  {
    double scale = fIntegral / hIntegral;
    h->Scale(scale);
  }

  h->Draw("hist");
  myfun.f->SetLineColor(kRed);
  myfun.f->Draw("SAME");
  myfun.c->Update();
  std::cout << "Deviazione standard: "
            << myfun.ssqm(h, myfun.f) << std::endl;
}

void rigenerazione(int M, int N, int nbins)
{
  static Function myfun;
  double xmin{0.};
  double xmax{4.};
  std::vector<std::vector<double>> all_counts(M, std::vector<double>(nbins, 0.));
  for (int j = 0; j < M; ++j)
  {
    TH1D *h = myfun.GH(N, nbins);
    for (int i = 1; i <= nbins; ++i)
    {
      all_counts[j][i - 1] = h->GetBinContent(i) / N;
    }
    delete h;
  }
  std::vector<double> mean(nbins, 0.);
  std::vector<double> sigma(nbins, 0.);

  for (int i = 0; i < nbins; ++i)
  {
    double sum = 0.;
    double sumsq = 0.;
    for (int j = 0; j < M; ++j)
    {
      sum += all_counts[j][i];
      sumsq += all_counts[j][i] * all_counts[j][i];
    }
    mean[i] = sum / M;
    double var = (sumsq - (M * mean[i] * mean[i])) / (M - 1);
    sigma[i] = std::sqrt(var);
  }
  std::vector<double> x(nbins, 0.), ex(nbins, 0.);
  for (int i = 0; i < nbins; ++i)
  {
    x[i] = xmin + (i + 0.5) * (xmax - xmin) / nbins;
  }

  TGraphErrors *g = new TGraphErrors(nbins, x.data(), mean.data(), ex.data(), sigma.data());
  g->SetTitle("incertezza di rigenerazione; x; conteggi medi");
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);
  g->SetLineColor(kBlue);

  double fintegral = myfun.f->Integral(xmin, xmax);
  double hintegral = std::accumulate(mean.begin(), mean.end(), 0.0);
  double scale = fintegral / hintegral;
  for (auto &m : mean)
    m *= scale;
  for (auto &s : sigma)
    s *= scale;

  static TCanvas *c2 = new TCanvas("c2", "Rigenerazione", 800, 600);
  g->Draw("AP");
  myfun.f->SetLineColor(kRed);
  myfun.f->Draw("SAME");
  c2->Update();
  double global_dev{0.};
  for (auto &s : sigma)
    global_dev += s;
  global_dev /= nbins;
  std::cout << "Deviazione standard media per bin: " << global_dev << std::endl;
}

int main(int argc, char **argv)
{
  TApplication app("app", &argc, argv);
  std::cout << "Avvio myMacro()..." << std::endl;

  prima_prova();
  rigenerazione(100, 1000000, 100);
  std::cout << "myMacro() terminata." << std::endl;
  app.Run();
  return 0;
}
