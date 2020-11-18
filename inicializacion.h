#ifndef __INICIALIZACION_H_
#define __INICIALIZACION_H_

#include <iostream>

using fptr = double(double, double, double, double);
using fptr2 = double(fptr, double, double, double);

const double XMIN = 0;
const double XMAX = 20;
const double dx = 0.5;
const double valor_inicial = 0.0;
//const double a_p = 5;
//const double b_p = 5;
const double eps = 0.001;

double caso_1(double x, double y, double a, double b);
double caso_2(double x, double y, double a, double b);

double euler(fptr fun, double x, double y, double h);
void integracion_euler(fptr fun, fptr2 alg, double h, double valor_inicial);

double rk4(fptr fun, double x, double y, double h);
void integracion_rk4(fptr fun, fptr2 alg, double h, double valor_inicial);

double cambio_max(fptr2 alg, fptr fun, double h, double valor_inicial);
double h_estable(fptr2 alg, fptr fun, double valor_inicial, double eps);

void max_global(fptr2 alg, fptr fun, double h, double valor_inicial);

double factorial(double n);

#endif
