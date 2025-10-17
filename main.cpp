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
class Function {
 public:
  double k{5.2};
  double phi{1.8};
  double b{0.2};

  TF1* f;
  TCanvas* c;

  Function() {
    f = new TF1("f", "pow(cos([0]*x + [1]),2) + [2]", 0., 4.);
    f->SetParNames("k", "phi", "b");
    f->SetParameters(k, phi, b);
    c = nullptr;
  }

  ~Function() {
    delete f;
    if (c) delete c;
  }

  TH1D* GH(int N, int nbins) {
    TH1D* h = new TH1D("h", "Distribuzione generata", nbins, 0., 4.);

    // Use TF1::GetRandom to draw x distributed according to f(x)
    for (int i = 0; i < N; ++i) {
      double x = f->GetRandom();
      h->Fill(x);
    }

    return h;
  }
  double ssqm(TH1D* h, TF1* f) {  // mi fa la deviazione standard
    int nbins = h->GetNbinsX();
    std::vector<double> sqm;
    sqm.reserve(nbins);
    for (int i = 0; i <= nbins; ++i) {
      sqm.push_back((h->GetBinContent(i) - f->Eval(h->GetBinCenter(i))) *
                    (h->GetBinContent(i) - f->Eval(h->GetBinCenter(i))));
    }

    double sum = std::accumulate(sqm.begin(), sqm.end(), 0.0);
    double dev = std::sqrt(sum / nbins);
    return dev;
  }

  void Draw() {
    c = new TCanvas("c", "Grafico cos^2(kx+phi)+b", 800, 600);
    f->SetTitle("f(x) = cos^{2}(k x + #varphi) + b; x; f(x)");
    f->Draw();
    c->Update();
  }
};

void prima_prova() {
  static Function myfun;  // serve per evitare che venga distrutta alla fine della funzione
  myfun.Draw();
  myfun.c->Update();
  TH1D* h = myfun.GH(1000000, 100);
  h->SetLineColor(kBlue);

  // Normalizza l'istogramma all'area della funzione
  double x_min = h->GetXaxis()->GetXmin();
  double x_max = h->GetXaxis()->GetXmax();
  double fIntegral = myfun.f->Integral(x_min, x_max);
  double hIntegral = h->Integral("width");
  if (hIntegral > 0.0 && fIntegral > 0.0) {
    double scale = fIntegral / hIntegral;
    h->Scale(scale);
  }

  h->Draw("hist");
  myfun.f->SetLineColor(kRed);
  myfun.f->Draw("SAME");
  myfun.c->Update();
  std::cout << "Deviazione standard quadratica media: "
            << myfun.ssqm(h, myfun.f) << std::endl;
}

int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  std::cout << "Avvio myMacro()..." << std::endl;

  prima_prova();

  std::cout << "myMacro() terminata." << std::endl;
  app.Run();
  return 0;
}
