#include "metod.hpp"
#include "funzioni.cpp"



int main(int argc, char **argv)
{
  try
  {
  TApplication app("app", &argc, argv);
  std::cout << "Avvio myMacro()..." << std::endl;
     gf::prima_prova();
  gf::rigenerazione(100, 1000000, 100);
  gf::bin_smearing(100, 100);
  gf::propagazione_parametri(100, 100);
  gf::fit_distribution();
  std::cout << "myMacro() terminata." << std::endl;
  app.Run();
  return 0;
  }
  catch (const std::invalid_argument& e) {
    std::cerr << "Argomento non valido: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
catch (const std::out_of_range& e) {
    std::cerr << "Out of range: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
catch (const std::runtime_error& e) {
    std::cerr << "Errore di runtime: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
catch (const std::exception& e) {
    std::cerr << "Eccezione std::exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
}
  
}
