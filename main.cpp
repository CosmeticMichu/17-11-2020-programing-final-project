#include <iostream>
#include "inicializacion.h"

int main(int argc, char **argv)
{
  double h1 = h_estable(rk4, caso_1, 0.0);
  integracion_rk4(caso_1, rk4, h1, 0.0);
  double h2 = h_estable(euler, caso_1, 0.0);
  integracion_euler(caso_1, euler, h2, 0.0);
  
  return 0;
}
