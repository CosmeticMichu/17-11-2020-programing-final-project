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
  int N = (XMAX - XMIN)/(k*h);
  for(int i = 0; i<=N; ++i)
  {
    double xi = XMIN + i*h*k;
    fout << xi << "\t" << aux << "\n";
    aux = euler(fun, xi, aux, k*h);
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
  int N = (XMAX - XMIN)/(k*h);
  for(int i=0; i<=N; ++i)
  {
    double xi = XMIN + k*h*i;
    fout << xi << "\t" << aux << "\n";
    aux = rk4(fun, xi, aux, k*h);
  }
  fout.close();
}

double cambio_max(fptr2 alg, fptr fun, double h, double valor_inicial)
{
  double aux = 0.0, aux1 = valor_inicial, aux2 = valor_inicial; // se usan tres auxiliares, el primero para guardar el cambio máximo, los otros dos para medirlo
  int N = (XMAX - XMIN)/(c*h);
  for(int i=0; i<=N; ++i)
  {
    double xi = XMIN + h*i*c;
    aux1 = alg(fun, xi, aux1, h*c); // calcula el siguiente paso de la funciòn
    if(std::fabs(aux2 - aux1)>aux)
    {
      aux = std::fabs(aux2 - aux1); //compara la diferencia entre los dos pasos y si es la mayor obtenida hasta el momento la almacena en el auxiliar
    }
    aux2 = aux1; 
  }
  return aux; // retorna el valor máximo
}

double h_estable(fptr2 alg, fptr fun, double valor_inicial, double eps)
{
  double h;
  for(int i=0; i<300; ++i)
  {
    h = std::pow(2, -i);
    if(cambio_max(alg, fun, h, valor_inicial)<=eps)
    {
      break;
    }
  }
return h;
}

void max_global(fptr2 alg, fptr fun, double h, double valor_inicial)
{
  int N = (XMAX - XMIN)/h;
  double aux = valor_inicial, aux1 = valor_inicial, xaux;
  for(int i = 0; i<=N; ++i)
  {
    double xi = XMIN + i*h;
    aux = alg(fun,xi, aux, h);
    if((aux - aux1)>0)
    {
      aux1 = aux;
      xaux = xi;
    }
  }
  std::cout << "El valor máximo de temperatura es " << aux1 << " y ocurre al tiempo " << xaux << "\n";
}
