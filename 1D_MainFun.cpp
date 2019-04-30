# include "1D_MainFun.h"

void Current_Power_calc(int,int,double *,double *,double *,double *,double *,double *,double *,double *,double *,double);
void VoltageSupply_correction(double,double,double *,double *,double *,char *,int dot);
void def_dtel_step(int,double,double *,double *,double,double,int,int,double *,double *,double *,double *);
void def_dtgas_step(int,double,double,double *,double,int,int,double *,double *,double *,double *,double *);
void mass_control(double *,double *,int,double *,double *);
double TRFaverage(double,double,double);

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Xi[Nmax][LEN+2],Pgas[LEN+2],Tgas[LEN+2],
           Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
           Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
           Er[LEN+1],Fir[LEN+2],PowExp;

    extern double tau,dt,WextL,WextR;//dte;
    extern int N,Nneg,Npos;
    extern bool TEEDF;

    int Nedf,Nchem;
    int i,n,nt,dot,dot1,dot2,dot3,dot4,Ndot1;
    double tic,etic,period = 0.0,Wext,ticRF = 0.0;
    double te_ddt = 1.0e-12,dte = te_ddt;

    init_read();
    init_gasDBparser();
    tic = init_data();
    etic = tic;

    Mesh_gen(Len);
    Mesh_GFcalc(Geom);

    Nedf = EEDF_read_CS(N,CSFile);
    Nchem = chem_make_react(Nedf,ChemFile);
    chem_read_react(Nchem,N);

    double Vext=3000.0/Eabs,dVext=0.0,PowCalc=0.0;
    double PowTav=0.0,IcalcTav=0.0;
    double Xi_Lav[N],Ni_Lav[N],Ng_Lav,Tg_Lav,Te_Lav,Tv_Lav,Jcond_Lav;
    //double Nel_Tav[LEN+2],Te_Tav[LEN+2],Jcond_Tav[LEN+2],Er_Tav[LEN+1];
    double Kch[LEN+2][Nchem],Sch[N+1][LEN+2],Wrad[LEN+2];//Kel[LEN+2][Nedf],
    double dTe[LEN+2],dTgas[LEN+2],dNel[LEN+2],EN[LEN+2],dEzN=0.0,dEz=0.0,Icalc=0.0;
    double Del[2][LEN+2],Muel[2][LEN+2],Di[N][LEN+2],Mui[N][LEN+2],Lam[LEN+2];
    double Eold[LEN+1],Jcond[LEN+2],Jdisp[LEN+2];
    double MNT,MNTr;

    Wext = fmax(WextL,WextR);

    mass_control(&MNT,&MNTr,N,&Ni[0][0],Tgas);
    gas_LenAverage(N,Geom,&Xi[0][0],&Ni[0][0],Ngas,Tgas,Te,Tv,Xi_Lav,Ni_Lav,&Ng_Lav,&Tg_Lav,&Te_Lav,&Tv_Lav);
    gas_LenPrint(&Xi[0][0],&Ni[0][0],N,Npos,Nneg,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Er,Fir,Jcond,tic);
    gas_TimePrint(Ni_Lav,N,Pgas[1],Tg_Lav,Ng_Lav,Te_Lav,Tv_Lav,0.0,Icalc,Vext,Wext,0.0,0.0,0.0,MNT,MNTr,tic,etic,dte,dt);

    for(i=1; i<=LEN; i++)
    {
        Wrad[i] = 0.0;
        Jcond[i] = 0.0;
    }

    double tg_chem,te_chem = 1e-8;
    double Ermax = 1.0e-10,Er0max = 1.0e-10;
    dot = 0, dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0, nt = 0;
    Ndot1 = 10;//int(Ndots/10);//150//20);///75);

    ///----------------------------------------------------------------------------------

    ///EEDF_Table-Data_Import/Calculation************************************************
    int err;
    err = EEDF_table_read(N,Nedf,EedfTblFl);
    if(err!=0)
        EEDF_table_calc(Nedf,N,LEN+2,&Ni[0][1],Ngas[1],Tgas[1],1.0e-5,EedfTblFl);

    ///----------------------------------------------------------------------------------

    ///*******************Main_Cycle_of_1D-RF-Discharge_Simulation*********************///
    do
    {
        dot++;
        dot1++;

        ///Ke_De_Mue_recalculation********************************************************
        for(i=1;i<=LEN;i++)
        {
            if(((nt==0) || (dot==Ndots) || (dVext>100) || (dot1==1)) && (Nel[i]>0.0))//|| || (dot1==Ndot1) || (dot1*dte/dt>0.005(dEzN>0.05)|| (dot1==1)
                EEDF_const_calc("E/N-case",Nedf,&Kch[i][0],Nel[i],Er[i],Ngas[i],&Te[i],&dTe[i],&Del[0][i],&Muel[0][i],tic);//&Kel[i][0]
                //EEDF_const_calc("Te-case",Nedf,&Kch[i][0],Nel[i],Er[i],Ngas[i],&Te[i],&dTe[i],&Del[0][i],&Muel[0][i],tic);//&Kel[i][0]

        }

        ///Transport_data_recalculation**************************************************
        if((nt==0) || (dot==Ndots)|| (dVext>100) || (dot1==1) || (dTgas[1]>10.0) || (dTe[1]*11605>100.0))///—равнение только в одной точке!!!???|| (dEz/Ez>0.001)
            Trasport_coefs_calc(N,Npos+Nneg,&Ni[0][0],Er,&Del[0][0],&Muel[0][0],&Di[0][0],&Mui[0][0],Lam,Pgas,Tgas,Ngas,Te);

        ///Chem_data_Coef-s_recalculation************************************************
        for(i=1;i<=LEN;i++)
        {
            if((nt==0) || (dot==Ndots) || (dTgas[i]>30.0) || (dTe[i]*11605>100.0))
                chem_const(&Kch[i][0],Nedf,Nchem,N,Te[i],Tgas[i],tic);
        }

        ///Heat-Transport_Eq_solve*******************************************************
        HeatTransport(Hgas,Ngas,Tgas,Lam,&Ni[0][0],&Xi[0][0],&Di[0][0],N,Er,Jcond,Wrad,dt);

        ///Poisson_Equation_solve********************************************************
        //Ermax = Poisson_SORsolve(Fir,Er,Vext,&Xi[0][0],&Ni[0][0],Ngas,Npos,Nneg,tic);
        Ermax = Poisson_SORsolve(Fir,Er,Vext,&Xi[0][0],&Ni[0][0],Ngas,Npos,Nneg,etic);

        ///Chem_Rates_recalculation******************************************************
        for(i=1;i<=LEN;i++)
            chem_runge_kutta4(i,&Ni[0][i],Nneg+Npos,N,&Kch[i][0],&Sch[0][i],Nchem,&Wrad[i],dt,dte,&tg_chem,&te_chem,tic,dot);

        ///Species-DDT+Chemistry_solve***************************************************
        ///+++Charged_species:
        for(n=0;n<=Npos+Nneg;n++)
            //SpecTransport(n,&Xi[n][0],&Ni[n][0],&Di[n][0],&Mui[n][0],Ngas,Er,&Sch[n][0],dte);
            SpecTransport_expl(&Xi[n][0],&Ni[n][0],&Di[n][0],&Mui[n][0],&Sch[n][0],Ngas,Er,dte);

        ///+++Neutral_species:
        for(n=Npos+Nneg+1;n<N;n++)
            //SpecTransport(n,&Xi[n][0],&Ni[n][0],&Di[n][0],&Mui[n][0],Ngas,Er,&Sch[n][0],dt);
            SpecTransport_expl(&Xi[n][0],&Ni[n][0],&Di[n][0],&Mui[n][0],&Sch[n][0],Ngas,Er,dt);

        ///Te-equation*******************************************************************
        TeTransport_expl(Te,Nel,&Del[1][0],&Muel[1][0],Jcond,Er,&Sch[N+1][0],dTe,dte);

        ///Gas-data_recalculation********************************************************
        for(i=1;i<=LEN;i++)
            gas_TP_calc(&Ni[0][i],&Xi[0][i],N,&Pgas[i],&Tgas[i],&dTgas[i],&Ngas[i],&Rogas[i],&Hgas[i],&Nel[i],&dNel[i],&Tv[i]);

        ///Boundary-conditions_recalculation*********************************************
        TransportBoundary_mod(N,Nneg+Npos,&Ni[0][0],&Xi[0][0],&Di[0][0],&Mui[0][0],Er,Ngas,Nel,Tgas,Pgas,Tv,Te,Tw);

        ///Power_recalculation***********************************************************
        Current_Power_calc(Npos,Nneg,&Mui[0][0],&Ni[0][0],Er,Eold,Fir,Jcond,Jdisp,&Icalc,&PowCalc,dt);

        tic += dt;
        etic += dte;
        ticRF += dte;

        ///RF-Period_averaging***********************************************************
        PowTav = TRFaverage(PowCalc,Wext,dte);
        IcalcTav = TRFaverage(Icalc,Wext,dte);

        ///Writing_and_saving_Gas-data***************************************************
        if(dot==Ndots)
        {
            mass_control(&MNT,&MNTr,N,&Ni[0][0],Tgas);

            gas_LenAverage(N,Geom,&Xi[0][0],&Ni[0][0],Ngas,Tgas,Te,Tv,Xi_Lav,Ni_Lav,&Ng_Lav,&Tg_Lav,&Te_Lav,&Tv_Lav);

            gas_LenPrint(&Xi[0][0],&Ni[0][0],N,Npos,Nneg,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Er,Fir,Jcond,tic);

            gas_TimePrint(Ni_Lav,N,Pgas[1],Tg_Lav,Ng_Lav,Te_Lav,Tv_Lav,Jcond_Lav,Icalc,Vext,Wext,PowCalc,PowTav,IcalcTav,MNT,MNTr,tic,etic,dte,dt);

            gas_SavePrint(&Ni[0][0],N,Pgas,Tgas,Nel,Te,Fir,tic);

            dot = 0;
        }

        ///Time-Step_Correction**********************************************************
        ///+++Charged_species:
        if((nt==0) || (dot3==1) || (fabs(1-Ermax/Er0max)>0.1))
        {
            def_dtel_step(dot3,tic,&dte,&te_ddt,te_chem,0.75,N,Npos+Nneg,Er,&Di[0][0],&Mui[0][0],&Ni[0][0]);//0.15

            Er0max = Ermax;

            dot3 = 0;
        }
        ///+++Neutral_species:
        if(nt==0 || dot2==Ndot1)
        {
            def_dtgas_step(dot2,tic,tg_chem,&dt,0.75,N,Npos+Nneg,&Xi[0][0],Rogas,Tgas,&Di[0][0],Lam);

            dot2 = 0;
        }

        ///External_Voltage-supply_correction********************************************
        if((nt==0) || (dot==Ndots) || (fabs(1.0-PowCalc/PowExp)>0.005) || (dVext>100))//&& (dot1==10)((fabs(1.0-PowTav/PowExp)>0.005) && (ticRF>1.0/Wext)))
        {
            VoltageSupply_correction(PowCalc,PowExp,&Vext,&dVext,Ngas,Geom,dot);//PowTav

            dot1 = 0;
        }

        nt++;
        dot2++;
        dot3++;

        if(ticRF>1.0/Wext)
        {
            PowTav = 0.0;
            IcalcTav = 0.0;
            ticRF = 0.0;
        }
    }
    while(tic<tau);

    return 0;
}
void Current_Power_calc(int Npos,int Nneg,double *Mui,double *Ni,double *E,double *Eold,double *Fi,double *Jcond,double *Jdisp,double *Icalc,double *PowCalc,double dt)
{
    int i,n;
    double Pow[LEN+2];

    for(i=1;i<=LEN;i++)
    {
        ///Conductivity_Current****************************
        Jcond[i] = Ni[i]*Mui[i];///electrons
        for(n=1; n<=Npos; n++) ///positive_ions
            Jcond[i] += Ni[n*(LEN+2)+i]*Mui[n*(LEN+2)+i];
        for(n=Npos+1; n<=Npos+Nneg; n++) ///negative_ions
            Jcond[i] += Ni[n*(LEN+2)+i]*Mui[n*(LEN+2)+i];
        Jcond[i] *= e*E[i];//[Abs/cm2*s]
        Jcond[i] = fabs(Jcond[i]);//[Abs/cm2*s]

        ///Displacement_Current****************************
        Jdisp[i] = (E[i]-Eold[i])/(dt*4.0*pi)*Eabs;

        ///Power_distribution******************************
        //Pow[i] = Jcond[i]*fabs(Fi[i]);
        Pow[i] = Jcond[i]*fabs(E[i]);///[abs/cm^3]
        //Pow[i] = (Jcond[i]+fabs(Jdis[i]))*fabs(Fi[i]);
    }

    ///Summarized_Current**********************************
    *Icalc = eKl/e*LENintegral(Jcond,Geom);//[Kl/s=A]

    ///Power_Supply_calculation****************************
    *PowCalc = eKl/e*LENintegral(Pow,Geom)*Eabs/Hght;///[W/cm^2]
}
void VoltageSupply_correction(double PowCalc,double PowExp,double *Vext,double *dVext,double *Ngas,char *Geom,int dot)
{
    int i;
    double dV,Vnew,
        Vmax=10000/Eabs;

    if(*dVext>100.0)
        *dVext = 0.0;

    double th = 1.05;
    double Vmin = 10.0;///[V]

    Vnew = *Vext;///[abs]
    if(fabs(1.0-(PowCalc/PowExp))>0.005)
    {
        if(PowCalc>PowExp)
            Vnew /= 1.005;
        else
            Vnew *= 1.005;
    }
    if(PowCalc/PowExp>th)
        Vnew /= th;

    ///Restriction_on_Vext*******************************************
    if(Vnew>Vmax)///[abs]
        Vnew = Vmax;

    if((Vnew!=*Vext) && dot==Ndots)
        printf("Vext = %.1lf[V]\n",Vnew);

    *dVext = fabs(Vnew - *Vext);
    *Vext = Vnew;

    /*
    double res;
    dEN = 0.0;
    if(!strcmp(Geom,"axial"))
    {
        for(i=1; i<=LEN; i++)
        {
            res = fabs(EN[i] - Enew/Ngas[i]*1e17*Eabs);//[Td]
            EN[i] = Enew/Ngas[i]*1e17*Eabs;//[Td]

            dEN += dEN*0.5*(l[i+1]+l[i])*(l[i+1]-l[i]);
        }
        dEN *= 2.0/(Len*Len);
    }
    else if(!strcmp(Geom,"cartesian"))
    {
        for(i=1; i<=LEN; i++)
        {
            res = fabs(EN[i] - Enew/Ngas[i]*1e17*Eabs);//[Td]
            EN[i] = Enew/Ngas[i]*1e17*Eabs;//[Td]

            dEN += dEN*(l[i+1]-l[i]);
        }
        dEN /= Len;
    }

    *dEzN += dEN;
    */
}
void def_dtel_step(int dot,double tic,double *dte,double *te_min0,double te_chem,double CFL,int N,int Nion,double *E,double *Di,double *Mui,double *Ni)
{
    double te_min, te_drift, te_diff=1.0, te_ddr=1.0;
    double dl,Cmax;

    Cmax = 0.15;//0.75;////0.75

    ///electronic_time_step***********************************************************
    te_drift = (l[1]-l[0])/fabs(0.5*(Mui[0]+Mui[1])*E[0]);//only_electrons
    for(int i=1; i<=LEN; i++)
    {
        dl = l[i+1]-l[i];

        te_diff = fmin(dl*dl/(2*Di[i]),te_diff);//only_electrons

        te_drift = fmin(dl/(fabs(0.5*(Mui[i+1]+Mui[i])*E[i])),te_drift);//only_electrons
        //te_drift = fmin(dl/(fabs(0.5*(Mui[i+1]+Mui[i])*E[i]*Ni[i]/Ni[0])),te_drift);//weighted

        te_ddr = fmin(fabs(dl/(2*Di[i]/dl-(0.5*(Mui[i+1]+Mui[i])*E[i]))),te_ddr);
    }

    te_min = fmin(te_diff,te_drift);
    te_min = fmin(te_ddr,te_min);
    te_min *= Cmax*CFL;
    //if(te_chem>=1.e-10)
        te_min = fmin(te_min,te_chem);

    ///dte_time-step_correction*******************************************************
    if(*dte<te_min)// && (Ermax/Er0max<0.99))//(fabs(dEr-dEr0)<0.1));//(dEr/dEr0<0.999))////(dEr/dEr0<0.95))//(fabs(dEr-dEr0)<0.1))//))(dEr<dEr0))
        *dte *= 1.01;
    if(te_min/(*te_min0)>=1.0001)// && (Ermax/Er0max<0.99))//(fabs(dEr-dEr0)<0.1));//(dEr/dEr0<0.999))////(dEr/dEr0<0.95))//(fabs(dEr-dEr0)<0.1))//))(dEr<dEr0))
        *dte *= 1.001;
    if(te_min/(*te_min0)<0.999 && *dte>te_min)
        *dte /= 1.005;//*te_min0/te_min;
    //if((1.e-12<dte<te_min) && (Ermax/Er0max>=1.01))//dEr/dEr0>=1.005))//(dEr-dEr0)>=0.05))//(dEr/dEr0>1.01))//(fabs(dEr-dEr0)>=0.1))//(dEr>dEr0))
    //dte *= 0.99;//1.005*dEr/dEr0;//*dEr/dEr0;//0.9975;
    //if((1.e-12<dte<te_min) && (Ermax/Er0max>=1.05))//dEr/dEr0>=1.005))//(dEr-dEr0)>=0.05))//(dEr/dEr0>1.01))//(fabs(dEr-dEr0)>=0.1))//(dEr>dEr0))
    //dte *= 0.99*Er0max/Ermax;//1.005*dEr/dEr0;//*dEr/dEr0;//0.9975;
    //if(dte>=0.1*te_min)
    //dte *= 0.999;
    if(*dte<1.e-12)
        *dte = 1.e-12;
    if(*dte>1.e-9)
        *dte = 1.e-9;
    if(*dte>te_min)
        *dte=te_min;

    *te_min0 = te_min;

    //*dte = 5.e-12;

    ///Logging************************************************************************
    FILE *Tmin;

    if(dot==0)
    {
        Tmin = fopen("Dte.txt", "w");
        fprintf(Tmin,"t,[s]\tdte\tte_min\tte_diff\tte_drift\tte_chem\tte_ddr\n");
    }
    else
        Tmin = fopen("Dte.txt", "a+");

    fprintf(Tmin,"%.10e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",tic,*dte,*te_min0,te_diff,te_drift,te_chem,te_ddr);
    fclose(Tmin);
}
void def_dtgas_step(int dot,double tic,double tg_chem,double *dtg,double CFL,int N,int Nion,double *Xi,double *Rogas,double *Tgas,double *Di,double *Lam)
{
    double tg_diff=1.0, tg_therm=1.0, tg_min;
    double Dmax,Cp,dl,Th;
    int i,ni,n;

    ///gaseous_time_step*************************************************************
    for(i=1; i<=LEN; i++)
    {
        dl = l[i+1]-l[i];
        gas_HCpSi_calc(Tgas[i],N);

        Cp = 0.0;
        Dmax = 0.0;
        for(n=1; n<N; n++)
        {
            ni = n*(LEN+2)+i;

            if(n>Nion)
                Dmax = fmax(Dmax,Di[ni]);

            Cp += HCpSi[1][n]*Xi[ni];//[эрг/г* ]
        }
        ///diffusional
        tg_diff = fmin(dl*dl/(2*Dmax),tg_diff);

        ///thermal
        Th = Lam[i]/(Rogas[i]*Cp);//Lam-[Ёрг/(с*см*K)]
        tg_therm = fmin(dl*dl/(2*Th),tg_therm);

    }
    tg_min = CFL*fmin(tg_therm, tg_diff);
    tg_min = fmin(tg_min,tg_chem);//tg_chem
    //tg_min *= CFL;

    ///dt_time-step_correction*******************************************************
    if(*dtg<tg_min)// && (Ermax/Er0max<0.99))//(fabs(dEr-dEr0)<0.1));//(dEr/dEr0<0.999))////(dEr/dEr0<0.95))//(fabs(dEr-dEr0)<0.1))//))(dEr<dEr0))
        *dtg *= 1.01;
    if(tg_min/(*dtg)<0.95)
        *dtg /= 1.05;//*te_min0/te_min;
    if(*dtg>=tg_min)
        *dtg=0.99*tg_min;
    if(*dtg<1.e-10)
        *dtg = 1.e-10;

    ///Logging************************************************************************
    FILE *Tmin;

    if(dot==0)
    {
        Tmin = fopen("Dt_gas.txt", "w");
        fprintf(Tmin,"t,[s]\tdt_gas\ttg_min\ttg_diff\ttg_therm\ttg_chem\n");
    }
    else
        Tmin = fopen("Dt_gas.txt", "a+");

    fprintf(Tmin,"%.7e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",tic,*dtg,tg_min,tg_diff,tg_therm,tg_chem);

    fclose(Tmin);
}
void mass_control(double *MNT,double *MNTr,int N,double *Ni,double *Tgas)
{
    int i,n,ni;
    double Sum,Sumr;

    Sum = 0.0;
    Sumr = 0.0;
    for(i=0; i<=LEN+1; i++)
    {
        for(n=1; n<N; n++)
        {
            ni = n*(LEN+2)+i;

            Sum += Mi[n]*Ni[ni];//*Tgas[i];
            Sumr += Mi[n]*Ni[ni]*(i+0.5);//*Tgas[i];
        }
    }
    Sumr /= (LEN+1)*(LEN+1);

    *MNT = Sum;
    *MNTr = Sumr;
}
bool converge(double *Y0, double *Y, double eps)
{
    double Norm, cnv;
    bool Conv;
    int i;

    cnv = 0.0;
    for(i=1; i<=LEN; i++)
    {
        //Norm = (Y[i]-Y0[i])*(Y[i]-Y0[i]);
        //Norm = sqrt(Norm);

        Norm = fabs(1-Y[i]/Y0[i]);

        if(Norm<eps)//(fabs(Fi[i]-Res)<exact);//1.0-Res/Fi[i])<1.0e-5
            cnv++;

        Y0[i] = Y[i];
    }

    if(cnv<LEN)
        Conv = false;
    else
        Conv = true;

    return Conv;
}
int sign(double x)
{
    int y;

    if(x>0)
        y = 1;
    else if(x<0)
        y = -1;
    else
        y = 0;

    return y;
}
double LENintegral(double *Y,char *Geom)
{
    int i;
    double Yy = 0.0;

    if(!strcmp(Geom,"axial"))
    {
        for(i=1; i<=LEN; i++)
            Yy += Y[i]*0.5*(l[i+1]+l[i])*(l[i+1]-l[i]);
        Yy *= 2*pi;//[—√—/с]
    }
    else if(!strcmp(Geom,"cartesian"))
    {
        for(int i=1; i<=LEN; i++)
            Yy += Y[i]*(l[i+1]-l[i]);
        Yy *= Hght;//[—√—/с]
    }

    return Yy;
}
double LENaverage(double *Y,char *Geom)
{
    double Yy;

    if(!strcmp(Geom,"axial"))
        Yy = LENintegral(Y,Geom)/(pi*Len*Len);
    else if(!strcmp(Geom,"cartesian"))
        Yy = LENintegral(Y,Geom)/(Len*Hght);

    return Yy;
}
double TRFaverage(double Y,double Wext,double dt)
{
    double Yav;

    Yav += Y*dt*Wext;

    return Yav;
}
