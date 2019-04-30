//Защита от повторного включения заголовка
#ifndef EEDF_CALC_H
#define EEDF_CALC_H
int EEDF_read_CS(int,char *);
void EEDF_calc(int,double *,double *,double *,double,double,double,double *,double *,double *,double *,double,double,int,double);//решение уравнения Больцмана
void EEDF_const_calc(char *,int,double *,double,double,double,double *,double *,double *,double *,double);
void EEDF_print(double *,double,double,double,double,char *);
void EEDF_table_calc(int,int,int,double *,double,double,double,char *);
void EEDF_table_print(int,int,char *);
int EEDF_table_read(int,int,char *);
#endif
