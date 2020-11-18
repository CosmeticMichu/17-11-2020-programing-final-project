#include <iostream>
#include "inicializacion.h"

int main(int argc, char **argv)
{
  /*
  double h1 = h_estable(rk4, caso_1, valor_inicial, eps);
  std::cout << h1 << "\n";
  integracion_rk4(caso_1, rk4, h1, valor_inicial);
  double h2 = h_estable(euler, caso_2, valor_inicial, eps);
  std::cout << h2 << "\n";
  integracion_euler(caso_2, euler, h2, valor_inicial);

max_global(rk4, caso_1, h1, valor_inicial);
max_global(euler, caso_2, h2, valor_inicial);
*/
  double h1 = h_estable(rk4, caso_1, valor_inicial, eps);
  integracion_rk4(caso_1, rk4, h1, valor_inicial);
  
  return 0;
}
