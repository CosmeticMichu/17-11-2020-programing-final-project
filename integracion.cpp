#include <iostream>
#include "inicializacion.h"
#include <fstream>
#include <cmath>

double caso_1(double x, double y, double a, double b)
{
  return a*std::pow(M_E, -2*x) - std::pow(y, 4); // el paràmetro b no se usa por lo que no afecta
}

double caso_2(double x, double y, double a, double b)
{
  return a*std::pow(M_E, -2*x) - b*(std::pow(y, 4) - 1);
}

double euler(fptr fun, double x, double y, double h)
{
return y + h*fun(x, y, a_p, b_p);
}

void integracion_euler(fptr fun, fptr2 alg, double h, double valor_inicial)
{
  std::ofstream fout("euler_method.txt");
  double aux = valor_inicial;
  int N = (XMAX - XMIN)/(h);
  for(int i = 0; i<=N; ++i)
  {
    double xi = XMIN + i*h;
    fout << xi << "\t" << aux << "\n";
    aux = euler(fun, xi, aux, h);
  }
  fout.close();
}

double rk4(fptr fun, double x, double y, double h)
{
  double k1, k2, k3, k4;
  k1 = fun(x, y, a_p, b_p);
  k2 = fun(x + 0.5*h, y + 0.5*k1*h, a_p, b_p);
  k3 = fun(x + 0.5*h, y + 0.5*k2*h, a_p, b_p);
  k4 = fun(x + h, y + h*k3, a_p, b_p);
  return y + h*(k1 + 2*k2 + 2*k3 + k4)/6;
}

void integracion_rk4(fptr fun, fptr2 alg, double h, double valor_inicial)
{
  std::ofstream fout("rk4.txt");
  double aux = valor_inicial;
  int N = (XMAX - XMIN)/(h);
  for(int i=0; i<=N; ++i)
  {
    double xi = XMIN + h*i;
    fout << xi << "\t" << aux << "\n";
    aux = rk4(fun, xi, aux, h);
  }
  fout.close();
}
