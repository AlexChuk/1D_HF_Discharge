# include "1D_MainFun.h"
# include "1DTransport.h"
# include "Gas_calc.h"
# include "1Dmesh.h"

extern double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];
double al_bound[Nmax+1][2],bet_bound[Nmax+1][2];

void Trasport_coefs_calc(int N,int Nion,double *Ni,double *E,double *Del,double *Muel,double *Di,double *Mui,double *Lam,double *Pgas,double *Tgas,double *Ngas,double *Te)
{
    double Damb,SumNi,DPT,MuPTi,EN;
    double A,B,xDiff,Dix,DN;
    int i,n,ni,k,ki,s,S;
    int Z[Nion+1];

    //Tgas[0] = Tgas[1];
    Te[0] = Te[1];
    Te[LEN+1] = Te[LEN];

    Z[0] = -1;
    for(n=1;n<=Nion;n++)
        Z[n] = sign(CMui[n]);

    for(i=0;i<=LEN+1;i++)
    {
        //Tgas[i] = 415.9;
        //Ngas[i] = 3.22e16*300/Tgas[i];
        EN = fabs(E[i])*1.0e17*Eabs/Ngas[i];///[Td]//5;//91.52;//

        for(n=0;n<N;n++)
        {
            ni = n*(LEN+2)+i;

            if(n==0)
            {
                Di[ni] = Del[i];
                Mui[ni] = Z[n]*Muel[i];//0.0;//

                *(Muel+ni+(LEN+2)) *= Z[n];///MuEel
                if(i==LEN+1)
                {
                    Di[0] = Di[1];
                    *(Del+(LEN+2)) = *(Del+(LEN+2)+1);
                    Di[LEN+1] = Di[LEN];
                    *(Del+2*(LEN+1)+1) = *(Del+2*(LEN+1));//*(Del+(LEN+2)+(LEN+1)) = *(Del+(LEN+2)+LEN);

                    Mui[0] = Mui[1];
                    *(Muel+(LEN+2)) = *(Muel+(LEN+2)+1);
                    Mui[LEN+1] = Mui[LEN];
                    *(Muel+2*(LEN+1)+1) = *(Muel+2*(LEN+1));
                }
            }
            else if(n>0 && n<=Nion)
            {
                DPT = pow(Tgas[i],CDi[n][1])/Pgas[i];
                MuPTi = 1.0;///

                ///Di(E/N)_calculation**********************************************
                Dix = CDi[n][0]*DPT;
                //DN = Dix*Ngas[i];
                if(CDi[n][3]!=0.0 && EN>=0.0)///Di(Ez/N)_dependence
                {
                    S = CDi[n][2]-1;

                    if(EN>CDi[n][S*2+3])
                        Dix  *= CDi[n][S*2+4];
                    else
                    {
                        s=0;
                        while(EN>CDi[n][s*2+3] && s<=S)
                            s++;

                        if(s==0)
                        {
                            A = (CDi[n][s*2+4]-1.0)/CDi[n][s*2+3];
                            B = 1.0;
                        }
                        else
                        {
                            A = (CDi[n][s*2+4]-CDi[n][(s-1)*2+4])/(CDi[n][s*2+3]-CDi[n][(s-1)*2+3]);
                            B = CDi[n][(s-1)*2+4]-A*CDi[n][(s-1)*2+3];
                        }
                        xDiff = A*EN+B;
                        Dix  *= xDiff;
                    }
                }
                Di[ni] = Dix;

                ///Mui(E/N)_calculation****************************************
                Mui[ni] = CMui[n]*MuPTi;
            }
            else
            {
                DPT = pow(Tgas[i],CDi[n][1])/Pgas[i];
                Di[ni] = CDi[n][0]*DPT;
                Mui[ni] = 0.0;//1.0;//e/me/Vm_av;
            }

            ///Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Эрг/(с*см*K)]Ref. - Nuttall_JRNBS_1957 (N2)
            if(n==N-1)
                //Lam[i] = 34*pow(Tgas[i],0.76);///[erg/cm*s*K]//[N2-Mankelevich_form_consists_with_Справочник-ФизВеличин]
                Lam[i] = 33*pow(Tgas[i],0.78);///[erg/cm*s*K]//[O2-My_approx_from_data_Справочник-ФизВеличин]
                //Lam[i] = 2.4/100*(6.025e-5*Tgas[i]+8.55e-3)*1e7;///[erg/cm*s*K]//[O2-Proshina_data]_=My*2.3
                //Lam[i] = 32*pow(Tgas[i],0.71);///[erg/cm*s*K]//[Ar-My_approx_from_data_Справочник-ФизВеличин]
        }
    }
}
void SpecTransport(int n,double *Xi,double *Ni,double *Di,double *Mui,double *Ngas,double *E,double *Si,double dt)
{
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    double Vdr[LEN+1];
    double *res;

    ///**************Defining_Drift-coefficient and Source_term******************
    int i;
    for(i=1;i<=LEN;i++)
    {
        if(i==1)
        {
            Vdr[i-1] = 0.5*(Mui[i]+Mui[i-1])*E[i-1];

            Si[i-1] = 0.0;
            Si[LEN+1] = 0.0;
        }

        Vdr[i] = 0.5*(Mui[i+1]+Mui[i])*E[i];

        //Si[i] = 0.0;
    }

    ///New_Jd=-D*N*dXi/dx_Vdr=N*Xi*********************************
    res = Transport_SWEEPsolve_mod(Ngas,Xi,Di,Vdr,Si,&al_bound[n][0],&bet_bound[n][0],dt);

    for(i=1;i<=LEN;i++)
    {
        Xi[i] = *(res+i);
        Ni[i] = *(res+i)*Ngas[i];///reNEW

        /*if(Ni[i]<1.0e-30)
            Ni[i] = 1.e-30;*/
    }
}
double* Transport_SWEEPsolve_mod(double *Ng,double *Xi,double *Di,double *Vi,double *Si,double *al_b,double *bet_b,double dt)
{
    /*
    //-----------------------------------------------------------------------------------------------------------------------------

	Transport equation (in Drift-Diffusion approximation ) solution with implicit SWEEP/Shuttle Method.

	dNi/dt = d(D*(dNi/dx))/dx - d(Vi*Ni)dx + S , S - source term.

	Drift-Diffusion Term (DDT) is implemented using exact solution for steady 1D DDT problem:

	d(RoV*F)/dx = dd(D*F)/dx*dx         (*),

	(F(x)-F(0))/(F(L)-F(0))=(exp(Pe*x/L)-1)/(exp(Pe)-1), Pe=RoV*L/D - Peclet Number

    In our Case:

    1)Species:
    J[i+1/2] = Vdr[i+0.5]*(N[i]-(N[i+1]-N[i])/(exp(Pe[i+0.5])-1));

    "Power Law Scheme" (Pantakar_1980_page-95) used for exp() approximation

    Linearization of Equation (*):

    J[i+1/2] = Vdr[i+0.5]*(N[i]-(1-0.1*Pe)^5(N[i+1]-N[i])/Pe), for 0<Pe<=5;

       [i-1]  [i]  [i+1]
    --|--x--|--x--|--x--|--

    A[i]*N[i] = A[i+1]*N[i+1] + A[i-1]*N[i-1]

    Coefficients (linear geometry):
    A[i+1] = D[i+0.5]*F(|Pe[i+0.5]|)+Max(-Vdr[i+0.5];0)
    A[i-1] = D[i-0.5]*F(|Pe[i-0.5]|)+Max(Vdr[i-0.5];0)
    A[i] = A[i+1]+A[i-1]+(Vdr[i+0.5]-Vdr[i-0.5])

    where:
    F(|Pe|) = Max(0;(1-0.1*|Pe|)^5)

    In "RADIAL" case:

    AA[i+1] = A[i+1]*r[i+0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i-1] = A[i-1]*r[i-0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i] = AA[i+1]+AA[i-1]+(r[i+0.5]*Vdr[i+0.5]-r[i-0.5]*Vdr[i-0.5])/(r[i]*(r[i+0.5]-r[i-0.5]))

    Chemistry part considered as explicit (source term)

	//-----------------------------------------------------------------------------------------------------------------------------
	*/

	double  D_L,D_R,V_L,V_R,F_L,F_R,dL,
            Pe,fPe,
            A_L,A_R,A_C;

    double  A,B,C,F,den,
            al[LEN+1],bet[LEN+1];

    ///SWEEP-SHUTTLE_CICLE*****************************************

    ///Defining_sweep_coefficients
    int i;
    for(i=1;i<=LEN;i++)
    {
        ///Defining_sweep_coefficients
        ///Diffusion
        D_L = D_R;
        if(i==1)
            D_L = 0.5*(Ng[i]+Ng[i-1])*0.5*(Di[i]+Di[i-1]);//0.5*(Di[i]*Ng[i]+Di[i-1]*Ng[i-1]);//
        D_R = 0.5*(Ng[i]+Ng[i+1])*0.5*(Di[i]+Di[i+1]);//0.5*(Di[i+1]*Ng[i+1]+Di[i]*Ng[i]);//

        ///Drift
        V_L = V_R;
        if(i==1)
            V_L = Vi[i-1]*0.5*(Ng[i]+Ng[i-1]);//0.5*(Vi[i]+Vi[i-1]);
        V_R = Vi[i]*0.5*(Ng[i+1]+Ng[i]);//0.5*(Vi[i+1]+Vi[i]);

        /*//Drift_NEW????
        V_L = V_R;
        if(i==1)
        {
            if(Vi[i-1]>0.0)
                V_L = Vi[i-1]*Ng[i-1];
            else
                V_L = Vi[i-1]*Ng[i];//0.5*(Vi[i]+Vi[i-1]);
        }
        if(Vi[i]>0.0)
            V_R = Vi[i]*Ng[i];
        else
            V_R = Vi[i]*Ng[i+1];
        */

        ///DDT_Coefficients:
        ///right-edge
        dL = 0.5*(l[i+2]-l[i]);
        F_R = V_R*dL;
        Pe = fabs(F_R/D_R);
        //Power-Law_scheme:
        fPe = pow((1.0-0.1*Pe),5.0);
        A_R = D_R*fmax(0.0,fPe)+fmax(-F_R,0.0);
        //Exp-scheme:
        //A_R = V_R/(exp(Pe)-1.0);

        ///left-edge
        dL = 0.5*(l[i+1]-l[i-1]);
        F_L = V_L*dL;
        Pe = fabs(F_L/D_L);
        //Power-Law_scheme:
        fPe = pow((1.0-0.1*Pe),5.0);
        A_L = D_L*fmax(0.0,fPe)+fmax(F_L,0.0);
        //Exp-scheme:
        //A_L = V_L/(exp(Pe)-1.0);

        ///Accounting for Geometry Factors_POWER-scheme:
        A_R = A_R*GF_R[i];
        A_L = A_L*GF_L[i];
        A_C = A_R+A_L+GF_R[i]*F_R-GF_L[i]*F_L;

        /*
        ///Accounting for Geometry Factors_EXP-scheme:
        A_R = A_R*GF_R[i];
        A_L = (A_L+V_L)*GF_L[i];
        A_C = (A_R+V_R)*GF_R[i]+A_L*GF_L[i];
        */

        ///SWEEP Coefficients_NEW-dt:
        A = -A_L*dt;///[i-1](left_cell)
        B = -A_R*dt;///[i+1](right_cell)
        C = A_C*dt+Ng[i];///[i](center_cell)
        F = Xi[i]*Ng[i]+Si[i]*dt;///RHS-part

        /*///SWEEP Coefficients_NEW-Stationary_do_not_work??:
        A = -A_L;///[i-1](left_cell)
        B = -A_R;///[i+1](right_cell)
        C = A_C;///[i](center_cell)
        F = Si[i];///RHS-part*/

        if(i==1)//see SpecTransport();
        {
            al[i-1] = al_b[0];///see_SpecTransportBoundary();
            bet[i-1] = bet_b[0];//0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;

        /*if(i==LEN)//see SpecTransport();
        {
            al[i] = al_b[1];///see_SpecTransportBoundary();
            bet[i] = bet_b[1];
        }*/

    }

    //Reverse_sweep_cycle************************************************
    for(i=LEN;i>=0;i--)
    {
        Xi[i] = al[i]*Xi[i+1]+bet[i];
        if(Xi[i]*Ng[i]<1.e-30)
            Xi[i] = 0.0;
    }

    return &Xi[0];
}
void TransportBoundary_mod(int N,int Nion,double *Ni,double *Xi,double *Di,double *Mui,double *E,double *Ngas,double *Nel,double *Tgas,double *Pgas,double *Tv,double *Te,double Twall)
{
    int n,nX=N-1,nR,nL;//N-1
    double Pe,Ped,D,Vt,Vd,JL,JR;
    double al,alXL,alXR,dlL,dlR,DXL,DXR,dNL,dNR;
    double Nnew,sign;

    dlL = 0.5*(l[2]-l[0]);
    dlR = 0.5*(l[LEN+2]-l[LEN]);

    ///Temperature_Boundary_conditions*****************************************
    ///Left_boundary**********************
    if(!strcmp(Geom,"axial"))
        Tgas[0] = Tgas[1];
    else if(!strcmp(Geom,"cartesian"))
        Tgas[0] = Twall;
    Ngas[0] = Pgas[0]/(kb*Tgas[0]);

    ///Right_boundary*********************
    Tgas[LEN+1] = Twall;
    Ngas[LEN+1] = Pgas[LEN+1]/(kb*Tgas[LEN+1]);

    ///Species_Boundary_conditions*********************************************
    dNL = 0.0;
    dNR = 0.0;
    for(n=0;n<N;n++)
    {
        nL = n*(LEN+2);
        nR = nL+LEN+1;

        ///Left_boundary***********************************************
        if(n<=Nion)///Charged_particles
        {
            if(Gamma[n][0] == 0.0)
            {
                Xi[nL] = Xi[nL+1];

                al_bound[n][0] = 1.0;
                bet_bound[n][0] = 0.0;
            }
            else if(Gamma[n][0] == 1.0)
            {
                Xi[nL] = 1.e-50;

                al_bound[n][0] = 0.0;
                bet_bound[n][0] = 0.0;
            }
            else
            {
                D = 0.5*(Di[nL+1]*Ngas[1]+Di[nL]*Ngas[0]);
                Vd = 0.5*(Mui[nL+1]*Ngas[1]+Mui[nL]*Ngas[0])*E[0];
                JL = -D*(Xi[nL+1]-Xi[nL])+Vd*Xi[nL+1];

                if(Xi[nL]-Xi[nL+1]>0.0 || JL>0.0)
                {
                    //Xi[nL] = Xi[nL+1];

                    //al_bound[n][0] = 1.0;
                }
                else
                {
                    Ped = dlR*fabs(Vd)/D;

                    if(n==0)
                        Vt = sqrt(8*kb*Te[0]*eV_K/(pi*Mi[n]));
                    else
                        Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]));
                    Vt  *= Ngas[0];
                    Pe = dlL*Vt/D;

                    al = (1.0+sign*0.25*Gamma[n][0]*Pe)/(1.0+Ped);

                    if(n==nX)
                        alXL = al;

                    Nnew = Xi[nL+1]/al;

                    if(n!=0)
                        dNL += 0.25*Gamma[n][0]*Vt*Xi[nL];//fabs(Nnew-Nold);//

                    Xi[nL] = Nnew;

                    al_bound[n][0] = 1.0/al;
                }

                bet_bound[n][0] = 0.0;
            }
        }
        else///Neutral_particles
        {
            if(Gamma[n][0] == 0.0)//Ni[n][0] = Ni[n][1];
            {
                Xi[nL] = Xi[nL+1];

                al_bound[n][0] = 1.0;
                bet_bound[n][0] = 0.0;

                if(n==nX)
                {
                    DXL = 0.5*(Di[nL+1]*Ngas[1]+Di[nL]*Ngas[0]);
                    alXL = 1.0;
                }
            }
            else///equation for DDT with kinetic wall flux
            {
                if(Xi[nL+1]-Xi[nL]>0.0)
                    sign = 1.0;
                else
                    sign = 0.0;

                D = 0.5*(Di[nL+1]*Ngas[1]+Di[nL]*Ngas[0]);
                if(n==nX)
                    DXL = D;

                Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]))*Ngas[0];
                Pe = dlL*Vt/D;

                al = (1.0+sign*0.25*Gamma[n][0]*Pe);

                if(n==nX)
                    alXL = al;

                Nnew = Xi[nL+1]/al;

                if(n!=nX)
                    dNL += 0.25*Gamma[n][0]*Vt*Xi[nL];//fabs(Nnew-Nold);//

                Xi[nL] = Nnew;

                al_bound[n][0] = 1.0/al;
                bet_bound[n][0] = 0.0;
            }
        }

        ///Right_boundary************************************************
        if(n<=Nion)///Charged_particles
        {
            if(Gamma[n][1] == 0.0)
            {
                Xi[nR] = Xi[nR-1];

                al_bound[n][1] = 1.0;
                bet_bound[n][1] = 0.0;
            }
            else if(Gamma[n][1] == 1.0)
            {
                Xi[nR] = 1.e-30;

                al_bound[n][1] = 0.0;
                bet_bound[n][1] = 0.0;
            }
            else
            {
                D = 0.5*(Di[nR]*Ngas[LEN+1]+Di[nR-1]*Ngas[LEN]);
                Vd = 0.5*(Mui[nR]*Ngas[LEN+1]+Mui[nR-1]*Ngas[LEN])*E[LEN];
                JR = -D*(Xi[nR]-Xi[nR-1])+Vd*Xi[nR-1];

                if(Xi[nR]-Xi[nR-1]>0.0 || JR<0.0)
                {
                    Xi[nR] = Xi[nR-1];

                    al_bound[n][1] = 1.0;
                }
                else
                {
                    Ped = dlR*Vd/D;

                    if(n==0)
                        Vt = sqrt(8*kb*Te[LEN+1]*eV_K/(pi*Mi[n]));
                    else
                        Vt = sqrt(8*kb*Tgas[LEN+1]/(pi*Mi[n]));
                    Vt  *=  Ngas[LEN+1];
                    Pe = dlR*Vt/D;

                    al = (1.0+0.25*Gamma[n][1]*Pe)/(1.0+Ped);

                    Nnew = Xi[nR-1]/al;

                    if(n!=0)
                        dNR += 0.25*Gamma[n][1]*Vt*Xi[nR];

                    Xi[nR] = Nnew;

                    al_bound[n][1] = al;
                }

                bet_bound[n][1] = 0.0;
            }
        }
        else///Neutral_particles
        {
            if(Gamma[n][1] == 0.0)//Ni[n][LEN+1] = Ni[n][LEN];
            {
                Xi[nR] = Xi[nR-1];

                al_bound[n][1] = 1.0;
                bet_bound[n][1] = 0.0;

                if(n==nX)
                {
                    DXR = 0.5*(Di[nR]*Ngas[LEN+1]+Di[nR-1]*Ngas[LEN]);

                    alXR = 1.0;
                }
            }
            else///equation for DDT with kinetic wall flux
            {
                if(Xi[nR]-Xi[nR-1]>0.0)
                    sign = 0.0;//-1.0;
                else
                    sign = 1.0;

                D = 0.5*(Di[nR]*Ngas[LEN+1]+Di[nR-1]*Ngas[LEN]);
                if(n==nX)
                    DXR = D;

                Vt = sqrt(8*kb*Tgas[LEN+1]/(pi*Mi[n]))*Ngas[LEN+1];
                Pe = dlR*Vt/D;

                al = (1.0+sign*0.25*Gamma[n][1]*Pe);

                if(n==nX)
                    alXR = al;

                Nnew = Xi[nR-1]/al;

                if(n!=nX)
                    dNR += 0.25*Gamma[n][1]*Vt*Xi[nR];

                Xi[nR] = Nnew;

                al_bound[n][1] = al;
                bet_bound[n][1] = 0.0;
            }
        }

        Ni[nL] = Xi[nL]*Ngas[0];
        Ni[nR] = Xi[nR]*Ngas[LEN+1];

        if(n==0)
        {
            Nel[nL] = Ni[nL];
            Nel[nR] = Ni[nR];
        }
    }

    ///for_Ground_state_n=nX:
    nL = nX*(LEN+2);
    nR = nL+LEN+1;
    if(dNL!=0.0)
    {
        ///OLD_version
        //bet_bound[n][0] = dNL*dlL/(Ngas[0]*Mi[n]*DXL*alXL);//

        ///NEW_version
        bet_bound[nX][0] = 0.0;//dNL*dlL/(DXL*alXL);//
        Xi[nL] = Xi[nL+1];//alXL+bet_bound[n][0];??????
        Ni[nL] = Xi[nL]*Ngas[0];
    }
    if(dNR!=0.0)
    {
        ///OLD_version
        //bet_bound[n][1] = -1.0*dNR*dlR/(Ngas[LEN+1]*Mi[n]*DXR);
        //Xi[n*(LEN+2)+LEN+1] += -1.0*bet_bound[n][1]/alXR;

        ///NEW_version
        bet_bound[nX][1] = -dNR*dlR/DXR;//
        Xi[nR] = (Xi[nR-1]-bet_bound[nX][1])/alXR;
        Ni[nR] = Xi[nR]*Ngas[LEN+1];
    }

    ///Ne-Te_Boundary_conditions***********************************************
    ///left_boundary**********************
    Te[0] = Twall/eV_K;//Te[1];
    Tv[0] = Tv[1];

    al_bound[N+1][0] = 1.0;
    bet_bound[N+1][0] = 0.0;

    ///Right_boundary**********************
    Te[LEN+1] = Twall/eV_K;//Te[LEN];
    Tv[LEN+1] = Tv[LEN];

    al_bound[N+1][1] = 1.0;
    bet_bound[N+1][1] = 0.0;

}
void HeatTransport(double *Hgas,double *Ngas,double *Tgas,double *Lam,double *Ni,double *Xi,double *Di,int N,double *E,double *J,double *Wrad,double dt)
{
    //Сетка по длине:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    //Учет нагрева электронами в упругих соударениях + нагрев газа полем + теплопроводность на стенку трубки
    //Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Эрг/(с*см^2*K)]Ref. - Nuttall_JRNBS_1957 (N2)
    //Hin += (Qel)*dt - Lam*(Tin-Tw)/pow(Rad/2.4,2.0)*dt;//w/o QE

    double  Lm,QL,QR,JHL,JHR,T0,Dn;
    int i,n,ni;

    //Calculation_of_enthalpy*************************************

    for(i=1;i<=LEN;i++)
    {
        //Учет теплопроводности******************

        QL = QR;
        if(i==1)
        {
            Lm = 0.5*(Lam[i]+Lam[i-1]);
            QL = - Lm*(Tgas[i]-Tgas[i-1]);
        }
        Lm = 0.5*(Lam[i+1]+Lam[i]);
        QR = - Lm*(Tgas[i+1]-Tgas[i]);

        //Учет диффузионного переноса************

        JHL = JHR;
        if(i==1)
        {
            T0 = (Tgas[i]+Tgas[i-1])*0.5;
            gas_HCpSi_calc(T0,N);

            JHL = 0.0;
            for(n=1;n<N;n++)//1
            {
                ni = n*(LEN+2)+i;

                Dn = 0.5*(Di[ni]*Ngas[i]+Di[ni-1]*Ngas[i-1]);
                JHL += - HCpSi[0][n]*Mi[n]*Dn*(Xi[ni]-Xi[ni-1]);
            }
        }
        JHR = 0.0;
        T0 = (Tgas[i+1]+Tgas[i])*0.5;
        gas_HCpSi_calc(T0,N);
        for(n=1;n<N;n++)//1
        {
            ni = n*(LEN+2)+i;

            Dn = 0.5*(Di[ni+1]*Ngas[i+1]+Di[ni]*Ngas[i]);
            JHR += - HCpSi[0][n]*Mi[n]*Dn*(Xi[ni+1]-Xi[ni]);
        }

        Hgas[i] += -dt*(GF_R[i]*(QR+JHR)-GF_L[i]*(QL+JHL))+dt*(fabs(J[i]*E[i])-Wrad[i]);///{J*E}=[erg/cm3/c]
    }
}
void SpecTransport_expl(double *Xi,double *Ni,double *Di,double *Mui,double *Sch,double *Ngas,double *E,double dt)
{
    int i;
    double Flux[LEN+1];
    double  D,Vdr,A,F,dL,Pe,fPe;

    ///Flux_calculations_Jd=-D*N*dXi/dx_Vdr=N*Xi*********************************
    for(i=1;i<=LEN;i++)
    {
        if(i==1)
        {
            ///DDT_Coefficients:
            //D = 0.5*(Di[1]*Ngas[1]+Di[0]*Ngas[0]);//0.5*(Ngas[1]+Ngas[0])*0.5*(Di[1]+Di[0]);//0.5*(Di[i]*Ng[i]+Di[i-1]*Ng[i-1]);//
            D = 0.5*(Ngas[1]+Ngas[0])*0.5*(Di[1]+Di[0]);

            //Vdr = 0.5*(Mui[1]*Ngas[1]+Mui[0]*Ngas[0])*E[0];
            Vdr = 0.5*(Mui[1]+Mui[0])*0.5*(Ngas[1]+Ngas[0])*E[0];

            /*dL = 0.5*(l[2]-l[0]);
            if(Vdr>=0.0)
                Flux[0] = -D*(Xi[1]-Xi[0])+Vdr*dL*Ngas[0]*Xi[0];
            else
                Flux[0] = -D*(Xi[1]-Xi[0])+Vdr*dL*Ngas[1]*Xi[1];
            */

            ///DDT_Flux:
            if(Vdr!=0.0)///DDT_Flux
            {
                dL = 0.5*(l[2]-l[0]);
                F = Vdr*dL;
                Pe = fabs(F/D);
                //Power-Law_scheme:
                fPe = pow((1.0-0.1*Pe),5.0);
                A = D*fmax(0.0,fPe)+fmax(F,0.0);

                Flux[0] = -A*(Xi[1]-Xi[0])+F*Xi[0];
            }
            else///Only_Diffusion_Flux
                Flux[0] = -D*(Xi[1]-Xi[0]);
        }

        ///DDT_Coefficients:
        //D = 0.5*(Di[i+1]*Ngas[i+1]+Di[i]*Ngas[i]);//0.5*(Ngas[i]+Ngas[i+1])*0.5*(Di[i]+Di[i+1]);
        D = 0.5*(Ngas[i]+Ngas[i+1])*0.5*(Di[i]+Di[i+1]);

        //Vdr = 0.5*(Mui[i+1]*Ngas[i+1]+Mui[i]*Ngas[i])*E[i];//0.5*(Mui[i+1]+Mui[i])*E[i];
        Vdr = 0.5*(Mui[i+1]+Mui[i])*0.5*(Ngas[i]+Ngas[i+1])*E[i];

        ///DDT_Flux:
        if(Vdr!=0.0)///DDT_Flux
        {
            dL = 0.5*(l[i+2]-l[i]);
            F = Vdr*dL;
            Pe = fabs(F/D);
            //Power-Law_scheme:
            fPe = pow((1.0-0.1*Pe),5.0);
            A = D*fmax(0.0,fPe)+fmax(-F,0.0);

            Flux[i] = -A*(Xi[i+1]-Xi[i])+F*Xi[i];
        }
        else///Only_Diffusion_Flux
            Flux[i] = -D*(Xi[i+1]-Xi[i]);
    }

    ///DDT_Equation_Solving****************************************************
    for(i=1;i<=LEN;i++)
    {
        Ni[i] += -dt*(GF_R[i]*Flux[i]-GF_L[i]*Flux[i-1])+dt*Sch[i];

        /*if(Ni[i]<1.0e-30)
            Ni[i] = 1.e-30;*/
    }
}
void TeTransport_expl(double *Te,double *Nel,double *DE,double *MuE,double *Je,double *E,double *SE,double *dTe,double dt)
{
    int i;
    double Flux[LEN+1];
    double Te0,D,Vdr,A,F,dL,Pe,fPe;
    double NE[LEN+2];

    ///NeTe_Recalculation********************************************************
    for(i=0;i<=LEN+1;i++)
        NE[i] = 3.0/2.0*Nel[i]*Te[i];///

    ///Flux_calculations_Jd=-D*N*dXi/dx_Vdr=N*Xi*********************************
    for(i=1;i<=LEN;i++)
    {
        if(i==1)
        {
            ///DDT_Coefficients:
            D = 0.5*(DE[1]+DE[0]);
            Vdr = 0.5*(MuE[1]+MuE[0])*E[0];///MuE<-->sign(-1)

            /*dL = 0.5*(l[2]-l[0]);
            if(Vdr>=0.0)
                Flux[0] = -D*(Xi[1]-Xi[0])+Vdr*dL*Ngas[0]*Xi[0];
            else
                Flux[0] = -D*(Xi[1]-Xi[0])+Vdr*dL*Ngas[1]*Xi[1];
            */

            ///DDT_Flux:
            if(Vdr!=0.0)///DDT_Flux
            {
                dL = 0.5*(l[2]-l[0]);
                F = Vdr*dL;
                Pe = fabs(F/D);
                //Power-Law_scheme:
                fPe = pow((1.0-0.1*Pe),5.0);
                A = D*fmax(0.0,fPe)+fmax(F,0.0);

                Flux[0] = -A*(NE[1]-NE[0])+F*NE[0];
            }
            else///Only_Diffusion_Flux
                Flux[0] = -D*(NE[1]-NE[0]);
        }

        ///DDT_Coefficients:
        D = 0.5*(DE[i]+DE[i+1]);
        Vdr = 0.5*(MuE[i+1]+MuE[i])*E[i];

        ///DDT_Flux:
        if(Vdr!=0.0)///DDT_Flux
        {
            dL = 0.5*(l[i+2]-l[i]);
            F = Vdr*dL;
            Pe = fabs(F/D);
            //Power-Law_scheme:
            fPe = pow((1.0-0.1*Pe),5.0);
            A = D*fmax(0.0,fPe)+fmax(-F,0.0);

            Flux[i] = -A*(NE[i+1]-NE[i])+F*NE[i];
        }
        else///Only_Diffusion_Flux
            Flux[i] = -D*(NE[i+1]-NE[i]);
    }

    ///DDT_Equation_Solving****************************************************
    for(i=1;i<=LEN;i++)
    {
        Te0 = Te[i];

        NE[i] += -dt*(GF_R[i]*Flux[i]-GF_L[i]*Flux[i-1]);//+dt*(fabs(Je[i]*E[i])/1.602e-12);//-SE[i]);///{J*E}=[erg/cm3/c]_1.602e-12;//[eV]=1.602e-12[erg]

        Te[i] = 2.0/3.0*NE[i]/Nel[i];///Te_Recalculation***********************

        dTe[i] = fabs(Te[i]-Te0);
    }
}
///******************************************NOT USED********************************************

double* Transport_SWEEPsolve(double *NNi,double *Di,double *Vi,double *Si,double *al_b,double *bet_b,double dt)
{
    /*
    //-----------------------------------------------------------------------------------------------------------------------------

	Transport equation (in Drift-Diffusion approximation ) solution with implicit SWEEP/Shuttle Method.

	dNi/dt = d(D*(dNi/dx))/dx - d(Vi*Ni)dx + S , S - source term.

	Drift-Diffusion Term (DDT) is implemented using exact solution for steady 1D DDT problem:

	d(RoV*F)/dx = dd(D*F)/dx*dx         (*),

	(F(x)-F(0))/(F(L)-F(0))=(exp(Pe*x/L)-1)/(exp(Pe)-1), Pe=RoV*L/D - Peclet Number

    In our Case:

    1)Species:
    J[i+1/2] = Vdr[i+0.5]*(N[i]-(N[i+1]-N[i])/(exp(Pe[i+0.5])-1));

    "Power Law Scheme" (Pantakar_1980) used for exp() approximation

    Linearization of Equation (*):

    J[i+1/2] = Vdr[i+0.5]*(N[i]-(1-0.1*Pe)^5(N[i+1]-N[i])/Pe), for 0<Pe<=5;

       [i-1]  [i]  [i+1]
    --|--x--|--x--|--x--|--

    A[i]*N[i] = A[i+1]*N[i+1] + A[i-1]*N[i-1]

    Coefficients (linear geometry):
    A[i+1] = D[i+0.5]*F(|Pe[i+0.5]|)+Max(-Vdr[i+0.5];0)
    A[i-1] = D[i-0.5]*F(|Pe[i-0.5]|)+Max(Vdr[i-0.5];0)
    A[i] = A[i+1]+A[i-1]+(Vdr[i+0.5]-Vdr[i-0.5])

    where:
    F(|Pe|) = Max(0;(1-0.1*|Pe|)^5)

    In "RADIAL" case:

    AA[i+1] = A[i+1]*r[i+0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i-1] = A[i-1]*r[i-0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i] = AA[i+1]+AA[i-1]+(r[i+0.5]*Vdr[i+0.5]-r[i-0.5]*Vdr[i-0.5])/(r[i]*(r[i+0.5]-r[i-0.5]))

    Chemistry part considered as explicit (source term)

	//-----------------------------------------------------------------------------------------------------------------------------
	*/

	double  D_L,D_R,V_L,V_R,dL,
            Pe,fPe,
            A_L,A_R,A_C;

    double  A,B,C,F,den,
            al[LEN+1],bet[LEN+1];

    //SWEEP-SHUTTLE_CICLE*****************************************

    //Defining_sweep_coefficients
    int i;
    for(i=1;i<=LEN;i++)
    {
        //Defining_sweep_coefficients
        //Diffusion
        D_L = D_R;
        if(i==1)
            D_L = 0.5*(Di[i]+Di[i-1]);
        D_R = 0.5*(Di[i+1]+Di[i]);

        //Drift
        V_L = V_R;
        if(i==1)
            V_L = Vi[i-1];//0.5*(Vi[i]+Vi[i-1]);
        V_R = Vi[i];//0.5*(Vi[i+1]+Vi[i]);

        //DDT_Coefficients:
        //right-edge
        dL = 0.5*(l[i+2]-l[i]);
        Pe = fabs(V_R*dL/D_R);
        fPe = pow((1.0-0.1*Pe),5.0);
        //A_R = D_R*std::max(0.0,fPe)+std::max(-V_R,0.0);
        A_R = D_R*fmax(0.0,fPe)+fmax(-V_R,0.0);

        //left-edge
        dL = 0.5*(l[i+1]-l[i-1]);
        Pe = fabs(V_L*dL/D_L);
        fPe = pow((1.0-0.1*Pe),5.0);
        //A_L = D_L*std::max(0.0,fPe)+std::max(V_L,0.0);
        A_L = D_L*fmax(0.0,fPe)+fmax(V_L,0.0);

        //Accounting for Geometry Factors:
        A_R = A_R*GF_R[i];
        A_L = A_L*GF_L[i];
        A_C = A_R+A_L+GF_R[i]*V_R-GF_L[i]*V_L;

        //SWEEP Coefficients:
        A = -A_L;///[i-1](left_cell)
        B = -A_R;///[i+1](right_cell)
        C = A_C+1.0/dt;///[i](center_cell)
        F = NNi[i]*1.0/dt+Si[i];///RHS-part

        if(i==1)//see SpecTransport();
        {
            al[i-1] = al_b[0];///see_SpecTransportBoundary();
            bet[i-1] = bet_b[0];//0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;

        if(i==LEN)//see SpecTransport();
        {
            al[i] = al_b[1];///see_SpecTransportBoundary();
            bet[i] = bet_b[1];
        }
    }

    //Reverse_sweep_cycle************************************************
    for(i=LEN;i>=0;i--)
    {
        NNi[i] = al[i]*NNi[i+1]+bet[i];
        if(NNi[i]<1.e-30)
            NNi[i] = 0.0;
    }

    return &NNi[0];
}
void TransportBoundary(int N,double *Ni,double *Xi,double *Di,double *Ngas,double *Tgas,double *Tv,double *Te,double Twall)
{
    int n,nX=11;
    double Pe,D,Vt,TsrL,TsrR;
    double al,alXL,alXR,dlL,dlR,DXL,DXR,dNL,dNR;
    double Nnew,sign,NgL,NgR;

    dlL = 0.5*(l[2]-l[0]);
    dlR = 0.5*(l[LEN+2]-l[LEN]);
    //TsrL = 0.5*(Tgas[0]+Tgas[1]);
    //TsrR = 0.5*(Tgas[LEN]+Tgas[LEN+1]);
    //NgL = 0.5*(Ngas[0]+Ngas[1]);
    //NgR = 0.5*(Ngas[LEN]+Ngas[LEN+1]);

    ///Boundary_conditions********************************************

    dNL = 0.0;
    dNR = 0.0;
    for(n=0;n<N;n++)
    {
       ///Left_boundary***********************************************
        if(Gamma[n][0] == 0.0)//Ni[n][0] = Ni[n][1];
        {
            Ni[n*(LEN+2)] = Ni[n*(LEN+2)+1];
            al_bound[n][0] = 1.0;
            bet_bound[n][0] = 0.0;

            if(n==nX)
            {
                DXL = 0.5*(Di[n*(LEN+2)+1]+Di[n*(LEN+2)]);
                alXL = 1.0;
            }
        }
        else///equation for DDT with kinetic wall flux
        {
            if(Ni[n*(LEN+2)+1]-Ni[n*(LEN+2)]>0.0)
                sign = 1.0;
            else
                sign = -1.0;

            D = 0.5*(Di[n*(LEN+2)+1]+Di[n*(LEN+2)]);
            if(n==nX)
                DXL = D;

            if(n==0)
                Vt = sqrt(8*kb*Te[0]*eV_K/(pi*Mi[n]));
            else
                Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]));

            //if(Vd!=0)
            //Vd = 0.5*(Mui[1]+Mui[0])*E;
            //Pe = Vd*0.5*(l[2]-l[0])/D;
            //al = (1.0+0.25*Gamma[n][0]*Pe*Vt/Vd)/(1.0+Pe);

            Pe = dlL*Vt/D;
            al = (1.0+sign*0.25*Gamma[n][0]*Pe);
            //al = (Tgas[0]/Tgas[1]+sign*0.25*Gamma[n][1]*Pe*TsrL/Tgas[1]);//Mankelevich_like_Jd=-D*N*dXi/dl

            if(n==nX)
                alXL = al;

            Nnew = Ni[n*(LEN+2)+1]/al;
            Ni[n*(LEN+2)] = Nnew;

            if(n!=nX)
                dNL += 0.25*Gamma[n][0]*Vt*Nnew*Mi[n];//fabs(Nnew-Nold);//

            al_bound[n][0] = 1.0/al;
            bet_bound[n][0] = 0.0;
        }

        ///Right_boundary************************************************
        if(Gamma[n][1] == 0.0)//Ni[n][LEN+1] = Ni[n][LEN];
        {
            Ni[n*(LEN+2)+LEN+1] = Ni[n*(LEN+2)+LEN];
            al_bound[n][1] = 1.0;
            bet_bound[n][1] = 0.0;

            if(n==nX)
            {
                DXR = 0.5*(Di[n*(LEN+2)+LEN+1]+Di[n*(LEN+2)+LEN]);
                alXR = 1.0;
            }
        }
        else///equation for DDT with kinetic wall flux
        {
            if(Ni[n*(LEN+2)+LEN+1]-Ni[n*(LEN+2)+LEN]>0.0)
                sign = -1.0;
            else
                sign = 1.0;

            D = 0.5*(Di[n*(LEN+2)+LEN+1]+Di[n*(LEN+2)+LEN]);
            if(n==nX)
                DXR = D;

            if(n==0)
                Vt = sqrt(8*kb*Te[LEN+1]*eV_K/(pi*Mi[n]));
            else
                Vt = sqrt(8*kb*Tgas[LEN+1]/(pi*Mi[n]));

            //if(Vd!=0)
            //Vd = 0.5*(Mui[LEN]+Mui[LEN+1])*E;
            //Pe = Vd*0.5*(l[LEN+2]-l[LEN])/D;
            //al = (1.0+0.25*Gamma[n][1]*Pe*Vt/Vd)/(1.0+Pe);//al_bound=Ni[LEN]/Ni[LEN+1]

            Pe = dlR*Vt/D;
            al = (1.0+sign*0.25*Gamma[n][1]*Pe);//al=Ni[LEN]/Ni[LEN+1]
            //al = (Tgas[LEN+1]/Tgas[LEN]+sign*0.25*Gamma[n][1]*Pe*TsrR/Tgas[LEN]);//Mankelevich_like_Jd=-D*N*dXi/dl

            if(n==nX)
                alXR = al;

            Nnew = Ni[n*(LEN+2)+LEN]/al;
            Ni[n*(LEN+2)+LEN+1] = Nnew;

            if(n!=nX)
                dNR += 0.25*Gamma[n][1]*Vt*Nnew*Mi[n];//fabs(Nnew-Nold);//

            al_bound[n][1] = al;
            bet_bound[n][1] = 0.0;
        }
    }

    ///for_N2[0]
    n = nX;
    if(dNL!=0.0)
    {
        ///Ni_version
        bet_bound[n][0] = dNL*dlL/(Mi[n]*DXL*alXL);
        Ni[n*(LEN+2)] += bet_bound[n][0];
    }

    if(dNR!=0.0)
    {
        ///Ni_version
        bet_bound[n][1] = -1.0*dNR*dlR/(Mi[n]*DXR);
        Ni[n*(LEN+2)+LEN+1] += -1.0*bet_bound[n][1]/alXR;
    }

    ///Temperature_Boundary_conditions*****************************************
    ///Left_boundary**********************
    if(!strcmp(Geom,"axial"))
        Tgas[0] = Tgas[1];
    else if(!strcmp(Geom,"cartesian"))
        Tgas[0] = Twall;

    ///Right_boundary*********************
    Tgas[LEN+1] = Twall;

    Ngas[0] = 0.0;
    Ngas[LEN+1] = 0.0;
    for(n=0;n<N;n++)
    {
        Ngas[0] += Ni[n*(LEN+2)];
        Ngas[LEN+1] += Ni[n*(LEN+2)+LEN+1];
    }

    for(n=0;n<N;n++)
    {
        Xi[n*(LEN+2)] = Ni[n*(LEN+2)]/Ngas[0];
        Xi[n*(LEN+2)+LEN+1] = Ni[n*(LEN+2)+LEN+1]/Ngas[LEN+1];
    }

    ///Ne-Te_Boundary_conditions***********************************************
    ///left_boundary**********************
    Te[0] = Te[1];
    Tv[0] = Tv[1];
    al_bound[N+1][0] = 1.0;
    bet_bound[N+1][0] = 0.0;

    ///Right_boundary**********************
    Te[LEN+1] = Te[LEN];
    Tv[LEN+1] = Tv[LEN];
    al_bound[N+1][1] = 1.0;
    bet_bound[N+1][1] = 0.0;

}
void TeTransport(double *Ne,double *Te,double *Lel,double *Vel,int Nal,double dt)
{
    double NeTe[LEN+2],De[LEN+2],Ve[LEN+1],Se[LEN+2];

    //Сетка по длине:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    NeTe[0] = 1.5*Ne[0]*Te[0];
    NeTe[LEN+1] = 1.5*Ne[LEN+1]*Te[LEN+1];

    //Defining_Drift-coefficient and Source_term******************
    int i;
    double dL;
    for(i=1;i<=LEN;i++)
    {
        NeTe[i] = 1.5*Ne[i]*Te[i];

        if(i==1)
        {
            De[0] = 2.0*Lel[0]/Ne[0]/3.0;
            De[LEN+1] = 2.0*Lel[LEN+1]/Ne[LEN+1]/3.0;
        }
        De[i] = 2.0*Lel[i]/Ne[i]/3.0;//0.666667

        dL = 0.5*(l[i+1]-l[i-1]);

        Ve[i-1] = 5.0*Vel[i-1]/3.0+0.5*(De[i]+De[i-1])*(Ne[i]-Ne[i-1])/dL;   //1.666667
        if(i==LEN)
        {
            dL = 0.5*(l[i+2]-l[i]);
            Ve[i] = 5.0*Vel[i]/3.0+0.5*(De[i]+De[i+1])*(Ne[i+1]-Ne[i])/dL;
        }

        Se[i] = 0.0;//(Jel[i]*E[i])-Ne*Sum(ki*Ni*dEei);
    }

    //проверить_геометрические факторы в поправке для сноса!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //Transport_equation_solve************************************
    double *res;

    res = Transport_SWEEPsolve(NeTe,De,Ve,Se,&al_bound[Nal][0],&bet_bound[Nal][0],dt);
    for(i=1;i<=LEN;i++)
        Te[i] = *(res+i)/(1.5*Ne[i]);
}
