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
  int N = (XMAX - XMIN)/h; // la cantidad de particiones del intervalo
  double aux=valor_inicial + 0.01; // un pequeño cambio en la condición inicial
  double aux1 = valor_inicial, aux2 = 0; 
  aux2 = std::fabs(alg(fun, 0, aux1, h) - alg(fun, 0, aux, h)); // guarda la diferencia del algoritmo implementado en el mismo punto con una pequeña alteración de la condición inicial 
  for(int i=0; i <= N; ++i)
  {
    double xi = 0 + i*h;
    aux = alg(fun, xi, aux, h);
    aux1 = alg(fun, xi, aux1, h); //estos guardan los puntos siguientes dados para cada uno de los casos (el de la condicion inicial normal y el de la condicion inicial un poco cambiada)
    if(std::fabs(alg(fun, xi, aux1, h) - alg(fun, xi, aux, h))<= aux2){
      aux2 = std::fabs(alg(fun, xi, aux1, h) - alg(fun, xi, aux, h)); // si en el siguiente punto la diferencia entre los dos es mayor, se remplaza el valor 
    }
    return aux2; // retorna la diferencia máxima en el intervalo 
  }
}

double h_estable(fptr2 alg, fptr fun, double valor_inicial, double eps)
{
  double h;
  for(int i=0; i<300; ++i)
  {
    h = std::pow(1.01, -i);
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
