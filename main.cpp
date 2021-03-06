#include <iostream>
#include "inicializacion.h"

int main(int argc, char **argv)
{
  
  //double h1 = h_estable(rk4, caso_2, valor_inicial, eps); // calcula el h para el que la diferencia es menor que un epsilon
 // std::cout << h1 << "\n";
  integracion_rk4(caso_2, rk4, 0.001, valor_inicial); // implementa el método usando el h; el k se usa para lo propuesto
  
  //double h2 = h_estable(euler, caso_1, valor_inicial, eps);
  //std::cout << h2 << "\n";
  
  integracion_euler(caso_1, euler, 0.001, valor_inicial);
  
  /*
  max_global(rk4, caso_2, h1, valor_inicial);
  */
  max_global(euler, caso_1, 0.001, valor_inicial);
  

  return 0;
}
