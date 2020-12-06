#include <iostream>
#include "inicializacion.h"
#include <fstream>
#include <cmath>

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
  for(int i=0; i<1000; ++i)
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
