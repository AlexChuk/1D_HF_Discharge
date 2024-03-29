# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
//# include <algorithm>

# include "Initials.h"
# include "1Dmesh.h"
# include "1DPoisson.h"
# include "1DTransport.h"
# include "Chemistry_calc.h"
# include "EEDF_calc.h"
# include "Gas_calc.h"

#ifndef MAINFUN_H
#define MAINFUN_H

# define LEN 200  //dots per axis
# define Nmax 20
# define CSmax 100
# define NRmax 200

const double
	pi = 3.141592653589,
	c = 2.997924562e+10,
	ma = 1.67e-24,//[�]
	me = 9.1e-28,//[�]
	e = 4.8e-10,//[���]
	eKl = 1.60217662e-19,//[��]
	Eabs = 300,//E[�/��]=300E[���.���]
	Na = 6.022e+23,
	kb = 1.38e-16,
	eV_K = 11605,//1 �� � ��������� �
	p0 = 1333.22, // ����-� �������� ���� --> ��� [���/��^3]
	exact = 1.0e-3,//�������� �����
    cm_eV = 8065.5447;//����-� �������� cm-1 --> ��


extern int NR,Ndots;//N,Nt,Nte,Nchem,Nedf,Ndots,
extern char Spec[Nmax][10],Spec_R[Nmax][10];
extern double Mi[Nmax],HCpSi[3][Nmax],CXi[Nmax][2][8];
extern double Emax,dE,dEev;
extern char Geom[10];
extern double Len,l[LEN+3],Hght;
extern double Gamma[Nmax][2],Tw;
extern double CDi[Nmax][25],CMui[Nmax];
extern char CSFile[50],ChemFile[50],EedfTblFl[50];
extern char RName[NRmax][100];

bool converge(double *,double *,double);
int sign(double);
double LENintegral(double *,char *);
double LENaverage(double *,char *);

#endif
