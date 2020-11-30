#ifndef __INICIALIZACION_H_
#define __INICIALIZACION_H_
//comentario de prueba
// Archivo de inicializacion del programa que incluye las constantes y la definición de las funciones

#include <iostream>

using fptr = double(double, double, double, double); // usado para llamar a los dos tipos de ecuaciones diferenciales
using fptr2 = double(fptr, double, double, double); // usado para llamar los algoritmos numéricos empleados

const double XMIN = 0; // valor del tiempo inicial
const double XMAX = 20; // valor del tiempo final
const double valor_inicial = 0.0; // valor de la temperatura en el tiempo inicial 
const double a_p = 5; // paràmetro de la ecuación diferencial 
const double b_p = 5; // paràmetro de la ecuación diferencial 
const double eps = 0.001; // epsilon usado para medir la estabilidad del mètodo 
const double k = 1; // paràmetro usado para no imprimir un número exagerado de datos en el .txt sino los necesarios
const double c = 1; // parámetro usado para optimizar la función h_estable

double caso_1(double x, double y, double a, double b); // definición de la ecuación diferencial de un parámetro, se introduce el parámetro b para poder usar también el fptr 
double caso_2(double x, double y, double a, double b); // definición de la ecuación diferencial de dos parámetros

double euler(fptr fun, double x, double y, double h); // función que hace un paso del método de euler
void integracion_euler(fptr fun, fptr2 alg, double h, double valor_inicial); // función que hace el paso a paso usando la función anterior y los imprime en un .txt 

double rk4(fptr fun, double x, double y, double h); // lo mismo que la anterior pero con rk4 
void integracion_rk4(fptr fun, fptr2 alg, double h, double valor_inicial);

double cambio_max(fptr2 alg, fptr fun, double h, double valor_inicial); // función que mide el cambio máximo del método  en todo el intervalo para un h especìfico
double h_estable(fptr2 alg, fptr fun, double valor_inicial, double eps); // función, que usa la función anterior, y retorna un h para el cual el cambio máximo es menor a un epsilon 

void max_global(fptr2 alg, fptr fun, double h, double valor_inicial); // funciòn que imprimer en la pantalla el màximo y el tiempo en el que ocurre


#endif
