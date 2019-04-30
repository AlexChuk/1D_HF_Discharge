//Защита от повторного включения заголовка
#ifndef CHEMISTRY_CALC_H
#define CHEMISTRY_CALC_H

int chem_make_react(int,char *);
void chem_read_react(int,int);
void chem_const(double *,int,int,int,double,double,double);
void chem_runge_kutta4(int,double *,int,int,double *,double *,int,double *,double,double,double *,double *,double,int);
void chem_half_implicit(double *,int,int,double *,int,double *,double,double *,double,int);
void chem_half_modimplicit(double *,int,int,double *,int,double *,double,double *,double,int);
double chem_sor_iterator(int,double,double,double,double,double,double,double);
double chem_newton_iterator(int,double,double,double,double,double,double);
void chem_spec_contrib(int,int,double,double);
int chem_VV_VT_make(int,char *,char *,char *,double *,int);
int chem_VV_VT_const(double *,int,int,char *,double,double);

#endif
