# include "1D_MainFun.h"
# include "Gas_calc.h"

double HCpSi[3][Nmax];

void gas_HCpSi_calc(double Tgas,int N)//нахождение внутренней энергии,энтальпии и теплоёмкостей по зад. температуре из аппроксиммаций Cp(T)
{
	double Hi,Cpi,Si,Mm,Cpt;
	int l,n,x;

	//T0 = Temp0;
	for(n=1;n<N;n++)
	{
		Mm = Mi[n];//[г]
		//Расчёт характеристик в точке T0
		if(Tgas>1000)
			x=0;
		else
			x=1;

		Cpi = 0;
		Hi = 0;
		Si = 0;
		for(l=0;l<5;l++)
		{
			Cpt = CXi[n][x][l]*pow(Tgas,l);

			Cpi += Cpt;//безразм
			Hi += Cpt/(l+1);//безразм

			if(l!=0)
				Si += Cpt/l;//безразм
		}
		Cpi = Cpi*kb/Mm;//[эрг/г*К]
		Hi = (Hi*kb*Tgas+CXi[n][x][5]*kb+CXi[n][x][7]*1.6e-12)/Mm;//[эрг/г]
		Si = (Si+CXi[n][x][0]*log(Tgas)+CXi[n][x][6])*kb/Mm;//[эрг/г*К]

		HCpSi[0][n] = Hi;//[эрг/г]
		HCpSi[1][n] = Cpi;//[эрг/г*К]
		HCpSi[2][n] = Si;//[эрг/г*К]
	}
}
void gas_TP_calc(double *Nn,double *Xn,int N,double *Pgas,double *Tgas,double *dTgas,double *Ngas,double *Rogas,double *Hgas,double *Nel,double *dNel,double *Tv)//расчёт температуры газа((H,P)-const)
{
	int n;
	double Xi[N],Nin,Roin,Tin,Pin,Hin,XMi,Xmi[N],hin;
	double Nout,Pout,Tout,Roout,Hout;
    double Hi,Cpi;
	double ftn,Ftn,Tnn,Tn,dT;

	*dNel = *Nel;

	//Присвоение значений из адресов переменных
	Tin = *Tgas;
	Pin = *Pgas;
	Hin = *Hgas;

	Nin = 0.0;//*Ngas;//0.0;
	for(n=0;n<N;n++)
		Nin += *(Nn+n*(LEN+2));//Nin += Ni[n*(LEN+2)];

	XMi = 0.0;
    for(n=0;n<N;n++)
    {
        Xi[n] = *(Nn+n*(LEN+2))/Nin;//*(Xn+n*(LEN+2));////*(Xn+n*(LEN+2));//*(Nn+n*(LEN+2))/Nin;//*(Xn+n*(LEN+2));//Xi[n] = Ni[n*(LEN+2)]/Nin;
        Xmi[n] = Xi[n]*Mi[n];
        XMi += Xmi[n];
    }
    Roin = XMi*Nin;

    //Расчёт температуры_методом Ньютона
	Tn = Tin;//300;
	//hin = Hin/Nin;
	Tnn = Tin;
	do
	{
		if(!Tnn==0)
            Tn = Tnn;

		Ftn = 0;
		ftn = 0;
		gas_HCpSi_calc(Tn,N);
		for(n=1;n<N;n++)//Npos+Nneg
		{
			Hi = HCpSi[0][n];//[эрг/г]
			Cpi = HCpSi[1][n];//[эрг/г*К]

			Ftn += Xmi[n]*Nin*Hi;//Roi1[n]*Hi;//Xi[n]*Hi;//
			ftn += Xmi[n]*Nin*Cpi;//Roi1[n]*Cpi;//Xi[n]*Cpi;//
		}

		//Ftn = Ftn*Nn - Hin;
		Ftn = Ftn - Hin;//hin;//

		//Tnn = Tn - Ftn/(ftn*Nn);
		Tnn = Tn - Ftn/ftn;

		dT = fabs(Tnn-Tn);

	}while(dT>0.001);//0.001
	Tout = Tnn;//Tin//300;
	if(Tout<300.0)
        Tout = 300.0;

	//Isobaric process**************************************
	Pout = Pin;

	//Return to using variables:
	Nout = Pout/(kb*Tout);
	//Nout = Nin;
    Roout = XMi*Nout;

	gas_HCpSi_calc(Tout,N);
	Hout = 0;
	for(n=1;n<N;n++)
		Hout += Xmi[n]*HCpSi[0][n];//[эрг/cm^3]
    Hout *= Nout;
    //Hout = Hin;//*Nout/Nin;

	//Return to using variables:
	*Pgas = Pout;
	*Tgas = Tout;
	*Ngas = Nout;
	*Rogas = Roout;
	*Hgas = Hout;
    *dTgas = fabs(Tin-Tout);

	for(n=0;n<N;n++)
    {
        *(Nn+n*(LEN+2)) = Nout*Xi[n];//Ni[n*(LEN+2)] = Nout*Xi[n];
        *(Xn+n*(LEN+2)) = Xi[n];
    }

    //концентрация электронов
    *Nel = Nout*Xi[0];
    *dNel = fabs(*dNel-*Nel);

    //Tv-calculation
    /*for(n=0;n<N;n++)
    {
        if(!strcmp(Spec[n],"N2(0)"))
            break;
    }*/
    n=11;
    *Tv = fabs(HCpSi[0][n+1]-HCpSi[0][n])*Mi[n]/kb/log((*(Nn+n*(LEN+2)))/(*(Nn+(n+1)*(LEN+2))));//*Tv = fabs(HCpSi[0][n+1]-HCpSi[0][n])*Mi[n]/kb/log(Ni[n*(LEN+2)]/Ni[(n+1)*(LEN+2)]);

	//************************************************************
}
void gas_LenPrint(double *Xi,double *Ni,int N,int Npos,int Nneg,double *Pgas,double *Tgas,double *Ngas,double *Hgas,double *Rogas,double *Nel,double *Te,double *Tv,double *Er,double *Fir,double *Jel,double tic)//запись в файл
{
	int i,k,n,nmax;
	FILE *log;
	int I[] = {1,int(LEN/2),LEN};
	double Ch;
	char ChL[][5] = {"_0","_L/2","_L"};

	nmax = N-1;//12;

	///Veusz_data_printing*****************************************************
	log = fopen("GUIData.txt", "w");

    fprintf(log,"r,cm\t");
    for(n=0;n<=nmax;n++)
        fprintf(log,"%s\t",Spec[n]);
    fprintf(log,"Charge\n");

    ///Concentrations
    for(i=0;i<LEN+2;i++)
    {
        //fprintf(log,"%.3lf\t",(i+0.5)*Len/(LEN+1));//old
        fprintf(log,"%.3lf\t",0.5*(l[i+1]+l[i]));
        for(n=0;n<=nmax;n++)
            fprintf(log,"%.2e\t",Ni[n*(LEN+2)+i]);

        Ch = -Ni[i] ;//electrons
        for(n=1;n<=Npos;n++)//positive ions
            Ch += Ni[n*(LEN+2)+i];
        for(n>Npos;n<=Npos+Nneg;n++)//negative ions
            Ch += -Ni[n*(LEN+2)+i] ;

        //fprintf(log,"%.2e\n",Ch);
        fprintf(log,"%.2e\n",fabs(Ch));
    }
    fprintf(log,"\n");

    ///Mole_fractions
    for(n=0;n<=nmax;n++)
        fprintf(log,"x(%s)\t",Spec[n]);
    fprintf(log,"\n");

    for(i=0;i<LEN+2;i++)
    {
        for(n=0;n<=nmax;n++)
            fprintf(log,"%.2e\t",Xi[n*(LEN+2)+i]);
        fprintf(log,"\n");
    }
    fprintf(log,"\n");

    //fprintf(log,"Tgas,K\tTv,K\tTe,K\tTe,eV\tE/N,Td\n");
    fprintf(log,"Ngas,cm-3\tTgas,K\tTe,K\tTe,eV\n");
    for(i=0;i<LEN+2;i++)
        //fprintf(log,"%.2e\t%.2lf\t%.1lf\t%.1lf\t%.3lf\t%.2lf\n",Ngas[i],Tgas[i],Tv[i],Te[i]*eV_K,Te[i],E*1.0e17*Eabs/Ngas[i]);
        fprintf(log,"%.2e\t%.2lf\t%.1lf\t%.3lf\n",Ngas[i],Tgas[i],Te[i]*eV_K,Te[i]);
    fprintf(log,"\n");

    fprintf(log,"Fi,V\tEr/N,Td\tEr,V/cm\t");
    for(i=0;i<LEN+1;i++)
        fprintf(log,"%.3e\t%.3e\t%.3e\n",Fir[i]*Eabs,Er[i]*Eabs*1.e17/Ngas[i],Er[i]*Eabs);
    fprintf(log,"%.3e\n",Fir[LEN+1]*Eabs);
    fprintf(log,"\n");

    /*//VDF_print
    fprintf(log,"n\t");
    for(i=0;i<3;i++)
        //fprintf(log,"VDF(%.2lfcm)\t",l[I[i]]);
        fprintf(log,"VDF%s\t",ChL[i],l[I[i]]);
    fprintf(log,"\n");
    for(n=11;n<N;n++)
    {
        fprintf(log,"%d\t",n-11);
        for(i=0;i<3;i++)
            fprintf(log,"%.2e\t",Ni[n*(LEN+2)+I[i]]/Ni[11*(LEN+2)+I[i]]);
        fprintf(log,"\n");
    }
    fprintf(log,"\n");*/

    fclose(log);
}
void gas_TimePrint(double *Ni,int N,double Pgas,double Tgas,double Ngas,double Te,double Tv,double Jel,double Icalc,double Vext,double Wext,double Pcalc,double PowTav,double IcalcTav,double MNT,double MNTr,double tic,double etic,double dte,double dt)//запись в файл
{
    int i,k,n;
    double NRF = 1.0;
	FILE *log;

    ///Output_on_monitor*******************************************************
	printf("Time = %.2e[s]\n",tic);
	printf("LENaveraged_Data:\n");
	printf("P = %.1lf[Torr]\tT = %.1lf[K]\tTv = %.1lf[K]\n",Pgas/p0,Tgas,Tv);
	printf("I = %.2lf[mA]\tVext = %.2lf[kV]\tPower = %.2e[W/cm^2]\n",Icalc*1e3,Vext*Eabs/1.e3,Pcalc);
	printf("Xe = %.2e\t%s = %.2e[cm-3]\t%s = %.2e[cm-3]\t%s = %.2e[cm-3]\n\n",Ni[0]/Ngas,Spec[0],Ni[0],Spec[1],Ni[1],Spec[4],Ni[4]);

    ///Output_to_file**********************************************************
	if(tic==0.0)
    {
        log = fopen("Gas_TimeData.txt", "w");

        //fprintf(log,"t,s\tdte,s\tdt,s\tPgas(t),Torr\tTgas(t),K\tTv(t),[K]\tTe(t),[K]\tEz(t),V/cm\tE/N(t),Td\tJel(t),[A/cm2]\tIel(t),[mA]\tNgas(t),[cm-3]\tXel\tMNT,[g*cm-3*K]\tMNTr,[g*cm-3*K]\t");
        fprintf(log,"t,s\tt_el,s\tdte,s\tdt,s\tNRF\tPgas(t),Torr\tTgas(t),K\tTe(t),[K]\tJel(t),[A/cm2]\tIel(t),[mA]\tVext(t),[kV]\tPower(t),[W/cm^2]\t<Iel>RF,[mA]\t<Power>RF,[W/cm^2]\tNgas(t),[cm-3]\tXel\tMNT,[g*cm-3*K]\tMNTr,[g*cm-3*K]\t");///w/o_Tv
        for(n=0;n<N;n++)//15
            fprintf(log,"%s(t)\t",Spec[n]);
        fprintf(log,"\n");
    }
    else
        log = fopen("Gas_TimeData.txt", "a+");

    if(Wext!=0)
    {
        Vext *= sin(2*pi*Wext*etic);
        NRF = etic*Wext;
    }

    //fprintf(log,"%.7e\t%.2e\t%.2e\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.2e\t%.1lf\t%.2e\t%.2e\t%.5e\t%.5e\t",tic,dte,dt,Pgas/p0,Tgas,Tv,Te*eV_K,E*Eabs,E*1.0e17*Eabs/Ngas,Jel*eKl/e,Icalc*1e3,Ngas,Nel/Ngas,MNT,MNTr);
    fprintf(log,"%.7e\t%.7e\t%.2e\t%.2e\t%.2lf\t%.1lf\t%.1lf\t%.1lf\t%.2e\t%.3lf\t%.3lf\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.5e\t%.5e\t",tic,etic,dte,dt,NRF,Pgas/p0,Tgas,Te*eV_K,Jel*eKl/e,Icalc*1.e3,Vext*Eabs/1.e3,Pcalc,IcalcTav*1.e3,PowTav,Ngas,Ni[0]/Ngas,MNT,MNTr);///w/o_Tv
    for(n=0;n<N;n++)//15
        fprintf(log,"%.2e\t",Ni[n]);
    fprintf(log,"\n");

	fclose(log);
}
void gas_SavePrint(double *Ni,int N,double *Pgas,double *Tgas,double *Nel,double *Te,double *Fir,double tic)//запись в save-файл
{
	int i,k,n;
	FILE *save;

    //запись параметров газа*******************************************************************
	save = fopen("Save_data.txt", "w");

	fprintf(save,"Time,sec\t%.2e\n",tic);

    //Data_print
    fprintf(save,"Pgas,Torr\t%.1lf\n",Pgas[i]/p0);

    fprintf(save,"Tgas,K\t");
    for(i=0;i<LEN+2;i++)
        fprintf(save,"%.1lf\t",Tgas[i]);
    fprintf(save,"\n");

    fprintf(save,"Te,eV\t");
    for(i=0;i<LEN+2;i++)
        fprintf(save,"%.2lf\t",Te[i]);
    fprintf(save,"\n");

    //Field
    fprintf(save,"Fi,V\t");
    for(i=0;i<LEN+2;i++)
        fprintf(save,"%.2e\t",Fir[i]);//(LEN+2) - due to 2 additional virtual points
    fprintf(save,"\n");

    //Ni[n]
    for(n=0;n<N;n++)
    {
        fprintf(save,"%s\t",Spec[n]);//,cm-3
        for(i=0;i<LEN+2;i++)
            fprintf(save,"%.2e\t",Ni[n*(LEN+2)+i]);//(LEN+2) - due to 2 additional virtual points
        fprintf(save,"\n");
    }

    //electron_data

    /*//VDF_print
    int I[] = {1,int(LEN/2),LEN};
    for(i=0;i<3;i++)
    {
        fprintf(save,"VDF(l=%.2lfcm)\t",l[I[i]]);
        for(n=11;n<N;n++)
            fprintf(save,"%.2e\t",Ni[n*(LEN+2)+I[i]]/Ni[11*(LEN+2)+I[i]]);
        fprintf(save,"\n");
    }
    fprintf(save,"\n");

    fprintf(save,"\n\n");*/

	fclose(save);
}
void gas_LenAverage(int N,char *Geom,double *Xi,double *Ni,double *Ngas,double *Tgas,double *Te,double *Tv,double *Xi_av,double *Ni_av,double *Ng_av,double *Tg_av,double *Te_av,double *Tv_av)
{
    int n;
    for(n=0;n<N;n++)
    {
        Xi_av[n] = LENaverage(&Xi[n*(LEN+2)],Geom);
        Ni_av[n] = LENaverage(&Ni[n*(LEN+2)],Geom);
    }

    *Ng_av = LENaverage(Ngas,Geom);
    *Tg_av = LENaverage(Tgas,Geom);
    *Te_av = LENaverage(Te,Geom);
    *Tv_av = LENaverage(Tv,Geom);
}


