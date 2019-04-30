#ifndef TRANSPORT_H
#define TRANSPORT_H

void Trasport_coefs_calc(int,int,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
double* Transport_SWEEPsolve(double *,double *,double *,double *,double *,double *,double);
double* Transport_SWEEPsolve_mod(double *,double *,double *,double *,double *,double *,double *,double);
void SpecTransport(int,double *,double *,double *,double *,double *,double *,double *,double);
void SpecTransport_expl(double *,double *,double *,double *,double *,double *,double *,double);
void HeatTransport(double *,double *,double *,double *,double *,double *,double *,int,double *,double *,double *,double);
void TeTransport(double *,double *,double *,double *,int,double);
void TeTransport_expl(double *,double *,double *,double *,double *,double *,double *,double *,double);
void TransportBoundary(int,double *,double *,double *,double *,double *,double *,double *,double);
void TransportBoundary_mod(int,int,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double);

#endif // TRANSPORT_H
