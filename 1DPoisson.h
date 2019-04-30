//Защита от повторного включения заголовка
#ifndef POISSON_CALC_H
#define POISSON_CALC_H

double Poisson_SORsolve(double *,double *,double,double *,double *,double *,int,int,double);
void Poisson_boundary(int,double *,int,double,double,double,double,double);

#endif
