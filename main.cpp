#include "metod.hpp"
#include "funzioni.cpp"




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
