#pragma once

#include"ode_solver.h"

double f(double x);
matrix fun(matrix x, matrix ud1, matrix ud2);
double MetodaFibonacciego(double a, double b, double epsilon);
double FiN(int n);
matrix df1(double t, matrix Y, matrix ud1, matrix ud2);
matrix fR(matrix x, matrix ud1, matrix ud2);

matrix fun2(matrix x0, matrix ud1, matrix ud2);
matrix df2(double t, matrix Y, matrix ud1, matrix ud2);
matrix fR2(matrix x, matrix ud1, matrix ud2);

matrix Fun3(matrix x1, matrix x2, matrix ud1);
bool g1(matrix x1);
bool g2(matrix x2);
bool g3(matrix x1, matrix x2, double alpha);
matrix df3(double t, matrix Y, matrix ud1, matrix ud2);
matrix fR3(matrix x, matrix ud1, matrix ud2);

matrix fun3(matrix x, matrix ud1, matrix ud2);

matrix fun4(matrix x, matrix ud1, matrix ud2);
matrix grad4(matrix x, matrix ud1, matrix ud2);
matrix hesj4(matrix x, matrix ud1, matrix ud2);

matrix fT4(matrix x, matrix ud1, matrix ud2);

matrix fR4(matrix x, matrix ud1, matrix ud2);
matrix gf(matrix x, matrix ud1, matrix ud2);

matrix fT5(matrix x, matrix ud1, matrix ud2);
matrix fR5(matrix x, matrix ud1, matrix ud2);

matrix f5(matrix x, matrix ud1, matrix ud2);
matrix f5_1(double a, matrix x, matrix ud1, matrix ud2);
matrix f5_2(double a, matrix x, matrix ud1, matrix ud2);

matrix fT6(matrix x, matrix ud1, matrix ud2);

matrix fR6(matrix x, matrix ud1, matrix ud2);
matrix df6(double t, matrix Y, matrix ud1, matrix ud2);