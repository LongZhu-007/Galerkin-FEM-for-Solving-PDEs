#ifndef COMMON_H
#define COMMON_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define eps 1e-4        /* 定义一个迭代求解的精度 */
#define itr_max 1000   /* 定义一个迭代求解的上限的步数 */
int Gmesh1D(double *meshes, int N, double x_min, double x_max);
double hat01(double x, double x0, double x1);
double hat02(double x, double x1, double x2);
double hat(double x, double x0, double x1, char d);
double D_hat01(double x, double x0, double x1);
double D_hat02(double x, double x1, double x2);
double D_hat(double x, double x0, double x1, char d);

double init_D_hats(double x0, double x1, double (*p)(double), char d1, char d2);
double init_hats(double x0, double x1, double (*q)(double), char d1, char d2);
double init_f_hat(double x0, double x1, double (*f)(double), char d);
int Stiff(double *KK, int N, double *meshes, double (*p)(double), double (*q)(double));
int right_b(double *F, int N, double *meshes, double (*f)(double));
int Boundary_element1(double *KK, double *meshes, int N, int *node, int m);
int Solver(double *KK, double *F, int N, int *nodes, int m, double *Bval, double *ui);
int UonNodes(double *res, double *meshes, int N, double *us, int *(nodes), int m, double *(val));
int Multipy(double *us, double *A, double *b, int Size);
double linfnorm(double *v, int Size);
void print_d_Array(double *Arr, int N);
void print_d_Matrix(double *Mat, int row, int col);

double p(double x);
double q(double x);
double f(double x);
#endif // COMMON_H
