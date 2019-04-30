# include "1D_MainFun.h"
# include "EEDF_calc.h"

# define NEmax 3000

int CStype = 10,CSout;
int CSref[Nmax][10][CSmax],CSR[CSmax][5];
double CS[NEmax][CSmax],Ith[CSmax];
double EedfTbl[300][CSmax+7];
double EN_Tbl,dEN_Tbl;
double Emax = 50.0;
double dEev = Emax/NEmax,
       dE = dEev*1.602e-12;//[eV]=1.602e-12[erg]

int EEDF_read_CS(int N,char *CSFile)//считывание сечений EEDF-процессов(возврат кол-ва реакций)
{
	FILE *cross;
    cross = fopen(CSFile, "r");

	FILE *log;
	log = fopen("Log_CS.txt", "w");
	fclose(log);

	log = fopen("Log_CS.txt", "a+");

	double X[150],Y[150],A[150],B[150];
	char CSstr[100],CSstr1[10],Cmt[100];
	double Stat,shift;

	char KeyW[][20] =
	{
		"ELASTIC",
		"ROTATION",
		"VIB-EXCITATION",
		"EXCITATION",
		"DISSOCIATION",
		"DEEXCITATION",
		"IONIZATION",
		"ATTACHMENT",
		"RECOMBINATION",
		"DISSATTACHMENT"
	};

	char symb[] = "------------------------------------------------------------";

	int i,k,K,j,J,err,n,n_t,Kw,add,d,s,Serr,rev;

	//считывание комментария
	fscanf(cross,"%s",&Cmt);
	do
	{
		fscanf(cross,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	j = 0;
	Serr = 0;
	fscanf(cross,"%s",&CSstr);//считывание первой строки
	while(strcmp(CSstr,"END")!=0)
	{
		add = 0;
		rev = 0;

		err = 0;
		for(i=0;i<CStype;i++)//сравнение с ключевыми словами
		{
			if(!strcmp(CSstr,KeyW[i]))
			{
				Kw = i;
				CSR[j][0] = Kw;//CS-type

				err = 1;
				break;
			}
		}
		if(err==0)//нет совпадений
		{
			printf("!!!Can't find Keyword!!! %s\n", CSstr);
			fprintf(log,"!!!Can't find Keyword!!! %s\n", CSstr);

			Serr += 1;
			break;
		}

		//считывание второй строки
		fscanf(cross,"%s",&CSstr);
		err = 0;
		for(n=0;n<N;n++)
		{
			if(!strcmp(CSstr,Spec[n]))
			{
				n_t = n;///target particle
				CSR[j][2] = n_t+1;
				err += 1;
				break;
			}
		}
		if(err==0)
		{
			printf("!!!Unknown target particle - %s in #%d cross section!!!\n ",CSstr,j);
			fprintf(log,"!!!Unknown target particle - %s in #%d cross section!!!\n ",CSstr,j);

			Serr += 1;
			///break;
		}

		///обработка 2-й строки
		if(Kw>=1)//w/o ELASTIC+ROT-EXCITATION
		{
			fscanf(cross,"%s%s",&Cmt,&CSstr);
			if(Kw==4 || Kw==9)
				fscanf(cross,"%s",&CSstr1);

			if(!strcmp(Cmt,"="))//forvard or rev
            {
                CSR[j][1] = -1;
				rev = 1;

				fscanf(cross,"%lf",&Stat);///Stat=g_gr/g_ex=st_w[nr_p]/st_w[nr_t];
            }
			if(!strcmp(Cmt,"->"))
                CSR[j][1] = 1;

			err = 0;
			for(n=0;n<N;n++)
			{
				if(!strcmp(CSstr,Spec[n]))
				{
					CSR[j][3] = n+1;//product particle
					err += 1;

					if(err==2)
						break;
				}
				if((!strcmp(CSstr1,Spec[n]))&&(Kw==4 || Kw==9))///"DISSOCIATION" OR "DISSATTACHMENT"
				{
					CSR[j][4] = n+1;//product particle
					err += 1;

					if(err==2)
						break;
				}
			}

			///предупреждение
			if((Kw==4 || Kw==9)&&(err<2))
			{
				printf("!!!Unknown products in CS: %s %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr,CSstr1);
				fprintf(log,"!!!Unknown products in CS: %s %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr,CSstr1);

				Serr += 1;
				///break;
			}
			else if((err==0)&&(Kw!=4 || Kw!=9))
			{
				printf("!!!Unknown product particle in CS: %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr);
				fprintf(log,"!!!Unknown product particle in CS: %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr);

				Serr += 1;
				///break;
			}
		}

		///обработка дополнения ко 2-й строке и считывание 3-ей строки
		fscanf(cross,"%s",CSstr);
		if(!strcmp(CSstr,"//"))//доп. реакции с таким же сечением "//add: N2(3) N2(5) //"
		{
			add = 0;
            fscanf(cross,"%s",CSstr);
			while(strcmp(CSstr,"//")!=0)
			{
				err = 0;
				for(n=0;n<N;n++)
				{
					if(!strcmp(CSstr,Spec[n]))
					{
						add += 1;
						CSR[j+add][0] = CSR[j][0];//type
						CSR[j+add][1] = CSR[j][1];//forvard or rev
						CSR[j+add][2] = n+1;//target particle
						CSR[j+add][3] = CSR[j][3];//product particle
						if(Kw==4)
                            CSR[j+add][4] = CSR[j][4];//product particle
						err += 1;
						break;
					}
				}

				if(err==0)
				{
					printf("!!!Unknown product particle in CS: %s -> add: %s - !!!\n ",Spec[n_t],CSstr);
					fprintf(log,"!!!Unknown product particle in CS: %s -> add: %s - !!!\n ",Spec[n_t],CSstr);

					Serr += 1;
					break;
				}

				fscanf(cross,"%s",CSstr);
			}

            ///сдвиг сечения или нет
            shift = 0;
			fscanf(cross,"%s",CSstr);
			if(!strcmp(CSstr,"shift"))
                shift = 1;

			fscanf(cross,"%lf%s",&Ith[j],Cmt);

			for(d=1;d<=add;d++)
			{
				int n0;
				n0 = CSR[j][2]-1;
				n = CSR[j+d][2]-1;
				if(shift>0)
                    Ith[j+d] = Ith[j] - (CXi[n][0][7]-CXi[n0][0][7]);///same product case!!
                else
                    Ith[j+d] = Ith[j];///без сдвига
			}
		}
		else
		{
			Ith[j] = atof(CSstr);
			if(Kw==1)//rotational
				Ith[j] = Ith[j]/11605;//[eV]

			fscanf(cross,"%s",Cmt);
		}

		///формирование массива ссылок на номер сечения
		for(d=0;d<=add;d++)
		{
			n_t = CSR[j+d][2]-1;

			s = CSref[n_t][Kw][0]+1;
			CSref[n_t][Kw][s] = j+d+1;//j+1??????????????????****************
			CSref[n_t][Kw][0] += 1;
		}

		///считывание 4-ей строки
		fscanf(cross,"%s",Cmt);

		///считывание и обработка сечения
		fscanf(cross,"%s",CSstr);
		if(!strcmp(CSstr,symb))
		{
			k = 0;
			//X[0] = 0.0;
			//Y[0] = 0.0;
			fscanf(cross,"%s",CSstr);
			while(strcmp(CSstr,symb)!=0)
			{
				X[k] = atof(CSstr);
				fscanf(cross,"%lf",&Y[k]);
				k += 1;
				fscanf(cross,"%s",CSstr);
			}
			K = k; //счётчик
			X[K] = 1500;//Emax;
			Y[K] = 0.0;//

			///Loging*****************************************************
			fprintf(log,"#CS=%d %s for particle %s\n",j,KeyW[Kw],Spec[n_t]);
			for(k=0;k<K;k++)
				fprintf(log,"%.2lf\t%.2e\n",X[k],Y[k]);
			fprintf(log,"\n");
			///Loging*****************************************************
		}
		else
		{
			printf("WARNING!Error while reading CS %s\t for particle %s\n",KeyW[Kw],Spec[n_t]);
			fprintf(log,"WARNING!Error while reading CS %s\t for particle %s\n",KeyW[Kw],Spec[n_t]);
			Serr += 1;
			break;
		}

		///линеаризация сечения
		for(k=0; k<K; k++)
		{
			A[k] = (Y[k+1]-Y[k])/(X[k+1]-X[k]);
			B[k] = Y[k]-A[k]*X[k];
		}
		//экстраполяция до лимита по энергии - сохранение наклона!!!!!!!!
		//A[K-1] = (A[K-2]+A[K-3]+A[K-4])/3;
		//B[K-1] = B[K-2];//(B[K-2]+B[K-3]+B[K-4])/3;

		///построение сечений по расчётной энергетической сетке
		k = 0;
		double Ei, CSj;
		for(i=0; i<NEmax; i++)
		{
			Ei = (i+0.5)*dEev;
			if(Ei<Ith[j])
				CSj = 0.0;
			else
			{
				if(Ei>X[k+1])//while(X[k+1]<Ei)
					k += 1;
				if(Ei>=X[K-1])
					k = K-1;

				CSj = A[k]*Ei+B[k];
				if(CSj<0)
				{
					//printf("WARNING!!!Cross section #%d is negative at point E=%.2lfeV\n",j,Ei);
					//fprintf(log,"WARNING!!!Cross section #%d is negative at point E=%.2lfeV\n",j,Ei);
					//Serr += 1;
					CSj = 0.0;
				}
			}

			CS[i][j] = CSj;///сечения,прочитанных из файла процессов, по точкам в центре ячеек сетки

			if(add!=0)
			{
                int I;
                for(d=1;d<=add;d++)
                {
                    if(Ei<Ith[j+d])
                        CS[i][j+d] = 0.0;
                    else
                    {
                        I = i - int((Ith[j]-Ith[j+d])/dEev);//сдвиг порога по энергии
                        if(I<0)
                        {
                            fprintf(log,"\nWARNING!!!Threshold shift is larger than Ith(target) - %s Cross section #%d \n",KeyW[Kw],j+d);
                            break;
                        }
                        CS[I][j+d] = CS[i][j];

                        if(i==NEmax-1)
                        {
                            while(I<=NEmax-1)
                            {
                                CS[I][j+d] = CS[i][j];//
                                I += 1;
                            }
                        }
                    }
                }

			}

		}

		///for_reverse_process***************************************
        if(rev!=0)
        {
            int nr_t = CSR[j][3]-1;
            int nr_p = CSR[j][2]-1;
            int Ii = int(Ith[j]/dEev-0.5);

            ///Stat=g_gr/g_ex=st_w[nr_p]/st_w[nr_t];
            for(i=0;i<NEmax;i++)
            {
                Ei = (i+0.5)*dEev;

                if(i+Ii<NEmax)
                    CS[i][j+1] = CS[i+Ii][j]*(1+Ith[j]/Ei)*Stat;
                else
                    CS[i][j+1] = CS[NEmax-Ii-1][j+1];
            }

            CSR[j+1][0] = 5;///type="DEEXCITATION"
            CSR[j+1][1] = 1;///
            CSR[j+1][2] = nr_t+1;///target
            CSR[j+1][3] = nr_p+1;///product
            Ith[j+1] = Ith[j];///threshold

            ///Поправка в массив ссылок
            s = CSref[nr_t][5][0]+1;
			CSref[nr_t][5][s] = j+2;
			CSref[nr_t][5][0] += 1;
        }

		///Loging*****************************************************
		/*fprintf(log,"CS#%d at calc-energy net\n",j);
		for(i=0; i<NEmax; i++)
		{
			fprintf(log,"%.2lf\t%.2e\n",dEev*(i+0.5),CS[i][j]);
			i += 50;
		}
		fprintf(log,"\n");
		fclose(log);
		*/
		///Loging*****************************************************

		//пока не понял зачем
		//CS[j][0] = 0.0;
		//if(X[1]<dE1) {CXC[as][k][0]=Y[1];}

		j += 1+add+rev;

		fscanf(cross,"%s",CSstr);
	}
	fclose(cross);

	//Printing diagnosic message*****************************************************
	printf("EEDF Cross Sections info:\n\n");
	if(Serr>0)
	{
		printf("ATTENTION!!!Number of errors = %d\n", Serr);
		fprintf(log,"ATTENTION!!!Number of errors = %d\n", Serr);
	}
	else
	{
		printf("Number of EEDF-processes = %d\n", j);
		fprintf(log,"Number of EEDF-processes = %d\n", j);
	}

	J = j;

	//Запись в Лог***********************************************************

	fprintf(log,"Log for CS-matrix\n\n");

	fprintf(log,"E,eV\t\t");
	for(j=0;j<J;j++)
		fprintf(log,"CSnum=%d\t\t",j+1);
	fprintf(log,"\n");

	for(k=0;k<NEmax;k++)
	{
		fprintf(log,"%.2lf\t",dEev*(k+0.5));
		for(j=0;j<J;j++)
			fprintf(log,"%.2e\t",CS[k][j]);
		fprintf(log,"\n");
		//k += 10;
	}
	fprintf(log,"\n\n");

	fprintf(log,"Log for CSR[][]-matrix\n\n");

	fprintf(log,"CSnum\t");
	for(i=0;i<5;i++)
		fprintf(log,"Ind=%d\t",i);
	fprintf(log,"\n\n");

	for(j=0;j<J;j++)
	{
		fprintf(log,"%d\t",j+1);
		for(i=0;i<5;i++)
			fprintf(log,"%d\t",CSR[j][i]);
		fprintf(log,"\n");
	}
	fprintf(log,"\n\n");

	int Nt = 12;
	fprintf(log,"Log for CSref[%s][][]-matrix\n\n",Spec[Nt]);

	fprintf(log,"CStype\t\t RSum for %s\t\t CS-num\n\n",Spec[Nt]);

	int jj;
	for(Kw=0;Kw<CStype;Kw++)
	{
		jj = CSref[Nt][Kw][0];//CSref[n_t][Kw][0];
		fprintf(log,"%s\t\t%d\t",KeyW[Kw],jj);
		if(jj>0)
		{
			for(j=1;j<=jj;j++)
				fprintf(log,"%d\t",CSref[Nt][Kw][j]);
		}
		fprintf(log,"\n");
	}
	fprintf(log,"\n\n");

	//**********добавление EEDF-реакций в файл кин.схемы**********************//

	int n0,n1,n2,n3;
	FILE *react;
	react = fopen("ReactionSet.txt", "w");

	fprintf(react,"%s\n",symb);
	fprintf(react,"\nEEDF_Reactions\n");
	fprintf(react,"\n%s\n",symb);
	jj = 0;
	for(j=0;j<J;j++)
	{
		n0 = CSR[j][0];//CS-type
		n1 = CSR[j][2]-1;//target
		n2 = CSR[j][3]-1;//product-1
		n3 = CSR[j][4]-1;//product-2

		if(n0==6)//ion
		{
			fprintf(react,"R%d\t e\t + %s\t ->\t  %s\t + e\t + e\t ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
			jj += 1;
		}
		else if(n0==7 || n0==8)//att+rec
		{
			fprintf(react,"R%d\t e\t + %s\t ->\t %s\t ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
			jj += 1;
		}
		else if(n0>1)
		{
			if(n3>0)
			{
                if(n0==9)//dissatt
                {
                    fprintf(react,"R%d\t e\t + %s\t ->\t %s\t + %s ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2],Spec[n3]);
                    jj += 1;
                }
                else
                {
                    fprintf(react,"R%d\t e\t + %s\t ->\t e\t + %s\t + %s ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2],Spec[n3]);
                    jj += 1;
                }
			}
			else if(n2>0 && n2!=n1)
			{
				fprintf(react,"R%d\t e\t + %s\t ->\t e\t + %s ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
				jj += 1;
			}
		}
	}
	//fprintf(react,"\n");
	//fprintf(react,"// %s\n",symb);

	fclose(react);

	CSout = J-jj;///Elast+Rott+V-to-Vexcit_are_excluded!!!

	printf("%d of %d EEDF-processes were added to ReactionSet.txt\n\n",jj,J);
	fprintf(log,"%d of %d EEDF-processes were added to ReactionSet.txt\n\n",jj,J);

	fclose(log);

	return jj;
}
void EEDF_calc(int N,double *Ne,double *Ni,double *Te,double E,double Tgas,double Nel,double *De,double *Mue,double *DE,double *MuE,double dte,double tic,int dot,double eps)//решение уравнения Больцмана
{
	int k,n,m,s,j,J,Jmax,nte;
    double Te0,Te1,Norm,E_kT;
    double Ee,Vdr,De_Muel;

	double A,B,C,F,den;
	double Vm[NEmax],Vmi[NEmax],CSNm[NEmax],D[NEmax],Vm_av,Del,Muel,DEel,MuEel;
	double Ur,Ul,Dr,Dl,Vr,Vl;

	//Сетка по энергии:
	/*
	 Ne(k)[0]	[1]   [2]	[3]		   		 [k-1]  [k]  [k+1]						  [NEmax-1]
		|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->E
	  k=0     1     2     3     4		   k-1    k    k+1                    			   NEmax
	*/

	nte = 0;
	Te1 = *Te;
	do
    //(nte=0;nte<20;nte++)
    {
        Te0 = Te1;

        //Difinition of D=A and Vm (Raizer)
        for(k=0;k<=NEmax-1;k++)
        {
            Vm[k] = 0.0;
            Vmi[k] = 0.0;
            CSNm[k] = 0.0;
            double V = 0.0;
            for(n=1;n<N;n++)//по всем компонентам смеси (с кот. сталк. эл-ны), кроме электронов
            {
                //elastic collisions
                Jmax = CSref[n][0][0];//m=0

                if(Jmax>0 && Ni[n]>0.0)
                {
                    for(j=1;j<=Jmax;j++)
                    {
                        J = CSref[n][0][j]-1;

                        //CS[k][J] = 1e-15;//Test_for_Druvestein
                        //V = 2.0e11;//Test_for_Maxwell

                        CSNm[k] += CS[k][J]*Ni[n];

                        V = CS[k][J]*sqrt(2.0/me*dE*(k+0.5))*Ni[n];//sqrt(2/me)=4.69e13

                        Vm[k] += V;//[1/s]
                        Vmi[k] += V/Mi[n];
                    }
                }
            }

            D[k] = 2*e*e*E*E/(3*me*Vm[k]);// 2*e*e/(3*me*E0*E0)=1875.46*E^2/Vm - [B/см]^2*с
            Vmi[k] = 2*me*Vmi[k];//Sum(2*me*Vi/Mi);
        }

        ///Цикл прогонки******************************************************
        ///определение прогоночных коэффициентов
        double al[NEmax],bet[NEmax];
        double Qinel,Qex;

        for(k=0;k<=NEmax-1;k++)
        {
            ///набор энергии от поля и упругие потери в соударениях
            if(k==0)
                Dl = D[k];
            else
                Dl = 0.5*(D[k]+D[k-1]);

            if(k==NEmax-1)
                Dr = D[k];
            else
                Dr = 0.5*(D[k+1]+D[k]);

            Vr = 0.5*(Vmi[k+1]+Vmi[k]);
            Vl = 0.5*(Vmi[k]+Vmi[k-1]);

            Ur = 0.5*Dr-Vr*(k+1)*dE;//point Ne[k]-right face
            Ul = 0.5*Dl-Vl*k*dE;//point Ne[k-1]-left face

            //without drift
            A = - Dl*k*dE/(dE*dE);//k-1(left)
            B = - Dr*(k+1)*dE/(dE*dE);//k+1(right)
            C = -(A + B) + 1/dte;//k(center)

            //with drift part
            if((Ur>=0)&&(Ul>=0))
            {
                A += -Ul/dE;
                C += Ur/dE;
            }
            if((Ur>=0)&&(Ul<0))
                C += (Ur-Ul)/dE;
            if((Ur<0)&&(Ul<0))
            {
                C += -Ul/dE;
                B += Ur/dE;
            }
            if((Ur<0)&&(Ul>=0))
            {
                A += -Ul/dE;
                B += Ur/dE;
            }

            ///inelastic collisions*************************************************************************************
            Qinel = 0.0;
            s = 0;
            for(n=1;n<N;n++)//суммирование по всем компонентам смеси (с кот. сталк. эл-ны)
            {
                if(Ni[n]>0.0)
                {
                    for(m=2;m<CStype;m++)//суммирование по всем типам процессов
                    {
                        Jmax = CSref[n][m][0];//число пар (Ntarget,CStype)

                        if(Jmax>0)//только по ненулевым типам процессов
                        {
                            if((m>=2)&&(m<=4))///excitation(rot,vib,elec,diss)
                            {
                                for(j=1;j<=Jmax;j++)
                                {
                                    J = CSref[n][m][j]-1;//номер сечения для пар (Ntarget,CStype)

                                    s = k+int(Ith[J]/dEev);
                                    if(s>=NEmax-1)
                                        s = NEmax-1;

                                    Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];

                                    Qinel += Qex*Ni[n];
                                }
                            }
                            if(m==5)///deexcitation
                            {
                                for(j=1;j<=Jmax;j++)
                                {
                                    J = CSref[n][m][j]-1;

                                    s = k-int(Ith[J]/dEev);

                                    Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k];
                                    if(s>0)
                                        Qex += CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];
                                    Qinel += Qex*Ni[n];
                                }
                            }
                            if(m==6)///ionization_2-models
                            {
                                for(j=1;j<=Jmax;j++)
                                {
                                    J = CSref[n][m][j]-1;

                                    /*//равнораспределение между 2-мя эл-нами********************************************
                                    s = 2*k+int(Ith[J]/dEev);
                                    if(s>=NEmax-1)
                                        s = NEmax-1;

                                    Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+2*CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];
                                    */
                                    //*********************************************************************************

                                    //ионизирующий эл-н уносит остаток, новый уносит E0-Райзер+Хагеллар****************
                                    E_kT = Tgas/eV_K;//E_kT=kT(Волошин)[эВ]//
                                    s = k+int((Ith[J]+E_kT)/dEev);
                                    if(s>=NEmax-1)
                                        s = NEmax-1;

                                    Qex = 0.0;
                                    if((k>=int(E_kT/dEev-0.5))&&((k<=int(E_kT/dEev+0.5))))
                                    {
                                        int kI = int((Ith[J]+E_kT)/dEev);
                                        int ki = 0;
                                        for(ki=kI;ki<NEmax;ki++)
                                            Qex += CS[ki][J]*sqrt(2/me*dE*(ki+0.5))*Ne[ki];
                                    }

                                    Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s]+Qex;
                                    //**********************************************************************************

                                    Qinel += Qex*Ni[n];
                                }
                            }
                            if(m>=7)///attachment+recomb+diss-attach
                            {
                                for(j=1;j<=Jmax;j++)
                                {
                                    J = CSref[n][m][j]-1;

                                    Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k];
                                    Qinel += Qex*Ni[n];
                                }
                            }


                            //rotation??
                        }
                    }
                }
            }
            ///*********************************************************************************************************

            ///RH-part
            F = Ne[k]/dte + Qinel;

            if(k==0)
            {
                al[k+1] = -B/C;
                bet[k+1] = F/C;
            }
            else if(k==NEmax-1)
            {}
            else
            {
                den = 1.0/(A*al[k]+C);
                al[k+1] = -B*den;
                bet[k+1] = (F-A*bet[k])*den;
            }
        }

        ///граничное условие**************************************************
        Ne[NEmax-1] = 0.0;
        //Ne[NEmax-1] = (F-A*bet[NEmax-1])/(A*al[NEmax-1]+C);

        ///цикл обратной прогонки*********************************************
        for(k=NEmax-2;k>=0;k--)
        {
            Ne[k] = al[k+1]*Ne[k+1]+bet[k+1];
            if(Ne[k]<1e-30)
                Ne[k] = 0.0;
        }

        ///Normalization***********************************************
        Norm = 0.0;
        for(k=0; k<NEmax; k++)
            Norm += Ne[k]*dEev;

        for(k=0; k<NEmax; k++)
            Ne[k] = Ne[k]*Nel/Norm;

        ///Ee_Te-calculation*******************************************
        Ee = 0.0;
        for(k=0;k<NEmax;k++)
            Ee += Ne[k]*(k+0.5)*dEev*dEev;
        Ee = Ee/Nel;//[eV]
        Te1 = 2.0/3.0*Ee;//[eV]

        nte++;

    }while(fabs(1.0-(Te1/Te0))>eps);// || (nte<100));//(fabs(Te1-Te0)>0.01);////(fabs(Te1-Te0)>0.001);

    ///Apply_to_used_matrices******************************************
    *Te = Te1;

    ///De_&_Mue-calculation********************************************
    Del = 0.0;
    for(k=0;k<NEmax;k++)
        Del += sqrt((k+0.5)*dE)*Ne[k]/CSNm[k];
    Del *= sqrt(2.0/me)/(3.0*Nel)*dEev;
    *De = Del;///[см2/с]

    Muel = 0.0;
    for(k=0;k<NEmax-1;k++)
        Muel += (k+1)*(Ne[k+1]/sqrt(k+1.5)-Ne[k]/sqrt(k+0.5))/CSNm[k];
    Muel *= -sqrt(2.0*dE/me)*e/(3.0*Nel)*dEev/dE;
    *Mue = Muel;///[cm2/c*abs]

    De_Muel = Del*Eabs/Muel;//[В]
    Vdr = Muel*E;///[см/с]
    //*Jel = e*Nel*Vdr;//[Abs/cm2*s]

    //Qel = me/Mi[11]*Ee*Vm_av*Nel*1.602e-12;//[erg/cm3*s]

    //QE = E*Jel*3.14*pi*Rad*Rad; //[abs*СГС/s]

    ///DE_&_MuE-calculation********************************************
    DEel = 0.0;
    for(k=0;k<NEmax;k++)
        DEel += pow((k+0.5)*dE,1.5)*Ne[k]/CSNm[k];
    DEel *= sqrt(2.0/me)/(3.0*Nel*Ee)*dEev*dEev/dE/Ee;
    *DE = DEel;///[см2/с]

    MuEel = 0.0;
    for(k=0;k<NEmax-1;k++)
        MuEel += (k+1)*(k+1)*(Ne[k+1]/sqrt(k+1.5)-Ne[k]/sqrt(k+0.5))/CSNm[k];
    MuEel *= -sqrt(2.0*dE/me)*e/(3.0*Nel*Ee)*dEev*dEev/dE/Ee;
    *MuE = MuEel;///[cm2/c*abs]


    ///Writing_data***************************************************
	//if(dot==Ndots)
		//EEDF_print(Ne,*Te,Nel,Norm,tic);

	//***************************************************************

}
void EEDF_const_calc(char *Txt,int Nedf,double *Kel,double Nel,double E,double Ngas,double *Te,double *dTe,double *De,double *Mue,double tic)//Kel(EEDF)-calculation
{
	int n,m,s,Smax,sL,sR;
    double Te0,Ki,Y[Nedf+5],dEN,A,B,EN,ddTe;

    Te0 = *Te;
    EN = fabs(E)/Ngas*1e17*Eabs;
    Smax = int(EN_Tbl/dEN_Tbl);

    ///From_EEDF_Table***********************************
    if(!strcmp(Txt,"E/N-case"))///E/N-case***************
    {
        if(EN>EN_Tbl)
        {
            for(m=0;m<Nedf+7;m++)
                Y[m] = EedfTbl[Smax][m+2];
        }
        else if(EN<EedfTbl[0][m+2])
        {
            for(m=0;m<Nedf+7;m++)
                Y[m] = EedfTbl[0][m+2];
        }
        else
        {
            s = int(EN/dEN_Tbl);
            if(EN>EedfTbl[s][0] && EN<EedfTbl[s+1][0])
                s++;

            dEN = dEN_Tbl;//EedfTbl[s][0]-EedfTbl[s-1][0];
            for(m=0;m<Nedf+7;m++)
            {
                A = (EedfTbl[s+1][m+2]-EedfTbl[s][m+2])/dEN;
                B = EedfTbl[s][m+2]-A*EedfTbl[s][0];
                Y[m] = A*EN+B;
            }

            /*
            if(s==0)
            {
                for(m=0;m<Nedf+7;m++)
                    Y[m] = EedfTbl[0][m+2];
            }
            else if(s>0 && EN<=EN_Tbl)
            {
                dEN = dEN_Tbl;//EedfTbl[s][0]-EedfTbl[s-1][0];
                for(m=0;m<Nedf+7;m++)
                {
                    A = (EedfTbl[s][m+2]-EedfTbl[s-1][m+2])/dEN;
                    B = EedfTbl[s-1][m+2]-A*EedfTbl[s-1][0];
                    Y[m] = A*EN+B;
                }
            }
            */
        }

        //*Te = Y[0];
        *dTe = fabs(*Te-Te0);
    }
    else if(!strcmp(Txt,"Te-case"))///Te-case**************!!!!!!!!CHECKKKKKKKKKKKKKKKKKKKKKKKKK!!!!!!!!!!!!!!!!!!
    {
        if(Te0>EedfTbl[Smax][2])
        {
            for(m=0;m<Nedf+7;m++)
                Y[m] = EedfTbl[Smax][m+2];
        }
        else if(Te0<EedfTbl[0][2])
        {
            for(m=0;m<Nedf+7;m++)
                Y[m] = EedfTbl[0][m+2];
        }
        else
        {
            sL = 0,sR = Smax;
            do
            {
                s = int((sL+sR)/2);

                if(Te0>EedfTbl[s][2])
                    sL = s;
                else
                    sR = s;
            }while(!((Te0>EedfTbl[s][2]) && (Te0<EedfTbl[s+1][2])));

            ddTe = EedfTbl[s+1][2]-EedfTbl[s][2];
            for(m=0;m<Nedf+7;m++)
            {
                A = (EedfTbl[s+1][m+2]-EedfTbl[s][m+2])/ddTe;
                B = EedfTbl[s][m+2]-A*EedfTbl[s][0];
                Y[m] = A*Te0+B;
            }
        }
    }

    ///Присвоение_значений используемым переменным
    *Mue = Y[1]/Ngas*Eabs;//Mue
    *De = Y[2]/Ngas;//1.5;//Proshina_check!
    *(Mue+(LEN+2)) = Y[3]/Ngas*Eabs;///MueE
    *(De+(LEN+2)) = Y[4]/Ngas;///DeE
    for(m=0;m<Nedf;m++)
        Kel[m] = Y[m+5];
    //Kel[16] /= 1.5;////Proshina_check!

    /*
    ///Logging*******************************************
    FILE *log;
    if(tic==0)
        log = fopen("Kel_data.txt", "w");
    else
        log = fopen("Kel_data.txt", "a+");

    fprintf(log,"t=%.2e[s]\n",tic);
	for(m=0; m<Nedf; m++)
	    fprintf(log,"%d\t%.2e\n",m+1,Kel[m]);
	fprintf(log,"\n");
	fclose(log);
	///Logging*******************************************
	*/
}
void EEDF_print(double *Ne,double Te,double Nel,double Norm,double tic)//запись EEDF в файл
{
	int k;

    /*
	printf("\nt = %.2e[s]\n",tic);
	printf("Te = %.2lf[eV]\t Ne = %.2e[cm-3]\n",Te,Nel);
	*/

	FILE *log;
	log = fopen("EEDF_data.txt", "w");

    fprintf(log,"t=%.2e[s]\t Te=%.2lf[eV]\n",tic,Te);
    fprintf(log,"Ne=%.2e[cm-3]\n",Nel);
	for(k=0; k<NEmax; k+=10)
        fprintf(log,"%.2lf\t%.2e\n",(k+0.5)*dEev,Ne[k]/Nel);
	fprintf(log,"\n\n");

	/*
	fprintf(log,"t=%.2e [s]\nNorm=%.2e [cm-3]\nTe=%.2lf [eV]\nE,[eV]\t Ne(E),[cm-3/eV]\n",dt*(nt+1),Norm,Te);//Nel
	for(k=0; k<NEmax; k+10)
		fprintf(log,"%.2lf\t%.2e\n",dEev*(k+0.5),Ne[k]/Nel);
	fprintf(log,"\n\n");
	*/
	fclose(log);
}
void EEDF_table_calc(int Nedf,int N,int NLen,double *Nn,double Ng,double Tg,double eps,char *EedfTblFl)
{
    double E,E_N,Te,Nel,Del,Muel,DEel,MuEel,E_Nmin,Vdr,Ki;
    double Kel[Nedf],Ni[N],Ne[NEmax];
    int st,m,n,k,I,dot = 0;

    E_Nmin = 0.5;///[Td]

    for(n=0;n<N;n++)
        Ni[n] = *(Nn+n*NLen);
    Nel = Ni[0];

    ///Initial_EEDF****************************************************
    Te = 0.5;
    for(k=0;k<NEmax;k++)
    {
        Ne[k] = 2*Nel*pow((k+0.5)*dEev/(Te*Te*Te*pi),0.5)*exp(-(k+0.5)*dEev/Te);//Maxwell
        if(Ne[k]<1e-30)
            Ne[k] = 0.0;
    }
	///****************************************************************

    EN_Tbl = 1000.0;
    dEN_Tbl = 5.0;
    E_N = E_Nmin;
    I = 0;
    st = 0;
    do
    {
        E = E_N*Ng*1e-17/Eabs;

        EEDF_calc(N,Ne,Ni,&Te,E,Tg,Nel,&Del,&Muel,&DEel,&MuEel,1.0e-11,0.0,dot,eps);

        for(m=0;m<Nedf;m++)
        {
            Ki = 0.0;
            for(k=0;k<NEmax;k++)//int(Ith[m])
                Ki += CS[k][m+CSout]*sqrt(2.0/me*dE*(k+0.5))*Ne[k]*dEev;
            Kel[m] = Ki/Nel;
        }

        //EEDF_print(Ne,Te,Nel,E_N);

        Vdr = Muel*E_N*Ng*1e-17/Eabs;

        ///EedfTable_Data_Saving*************************************************
        EedfTbl[I][0] = E_N;
        EedfTbl[I][1] = Vdr;
        EedfTbl[I][2] = Te;
        EedfTbl[I][3] = Muel*Ng/Eabs;//[1/(cm*s*V)]
        EedfTbl[I][4] = Del*Ng;//[1/(cm*s)]
        EedfTbl[I][5] = MuEel*Ng/Eabs;//[1/(cm*s*V)]
        EedfTbl[I][6] = DEel*Ng;//[1/(cm*s)]
        for(int m=0;m<Nedf;m++)
            EedfTbl[I][m+7] = Kel[m];

        if(I==0)
            E_N = 0.0;
        E_N += dEN_Tbl;

        if(E_N>EN_Tbl && st!=0)
        {
            E_N = EN_Tbl;
            st++;
        }

        I++;

    }while(E_N<=EN_Tbl);

    EEDF_table_print(Nedf,I,EedfTblFl);
}
void EEDF_table_print(int Nedf,int I,char *EedfTblFl)
{
    FILE *tbl;
    int i,m;

    tbl = fopen(EedfTblFl, "w");
    fprintf(tbl,"E/N,Td\tVdr,cm/s\tTe,eV\tMue*N,1/(cm*s*V)\tDe*N,1/(cm*s)\t\tMuE*N,1/(cm*s*V)\tDE*N,1/(cm*s)\t");
    for(m=0;m<Nedf;m++)
        fprintf(tbl,"%s\t",RName[m]);
    fprintf(tbl,"\n");

    for(i=0;i<I;i++)
    {
        fprintf(tbl,"%.2lf\t",EedfTbl[i][0]);
        for(m=1;m<Nedf+7;m++)
            fprintf(tbl,"%.3e\t",EedfTbl[i][m]);
        fprintf(tbl,"\n");
    }

    printf("EEDF_Table_has_been_created_File:O2_EEDF_Ke_Table.txt\n\n");

    fclose(tbl);
}
int EEDF_table_read(int N,int Nedf,char *EedfTableFile)
{
    FILE *tbl;
    tbl = fopen(EedfTableFile, "r");

    FILE *log;
	log = fopen("Log_EedfTable.txt", "w");

    char Rstr[50];
    int i,n,err;

    if(!tbl)
    {
        err++;

        printf("\nNO EEDF-data file with %s name!!\n",EedfTableFile);
        printf("EEDF-data will be recalculated!\n");

        fprintf(log,"NO EEDF-data file with %s name!!\n",EedfTableFile);
        fprintf(log,"EEDF-data will be recalculated!\n");
    }
    else
    {
        fprintf(log,"The EEDF-file %s has been scanned\n",EedfTableFile);
        ///Scanning until the name of processes********************************************
        do
        {
            fscanf(tbl,"%s",Rstr);
        }while(strcmp("DE*N,1/(cm*s)",Rstr));

        ///Scanning names of processes and check them on matching with CS-data*************
        err = 0;
        for(n=0;n<Nedf;n++)
        {
            fscanf(tbl,"%s",Rstr);
            if(strcmp(Rstr,RName[n]))
            {
                fprintf(log,"Reaction R#%d - %s is not matching the e-Process %s from CS-file\n",n+1,Rstr,RName[n]);
                err++;
            }
        }

        ///Saving_the_data_from_EEDF-Table*************************************************
        if(err==0)
        {
            i = 0;
            do
            {
                for(n=0;n<Nedf+7;n++)
                    fscanf(tbl,"%lf",&EedfTbl[i][n]);

                EN_Tbl = fmax(EN_Tbl,EedfTbl[i][0]);

                i++;
            }while(i<300);
            dEN_Tbl = EedfTbl[2][0] - EedfTbl[1][0];

            printf("\nThe EEDF-data has been Successfully imported from %s!!!\n\n",EedfTableFile);
            fprintf(log,"The EEDF-data has been Successfully imported!!!\n");
        }
        else
        {
            printf("\nWarning!!!Problems while EEDF-data import have been detected!See Log_EedfTable.txt\n");
            printf("EEDF-data will be recalculated!\n\n");
            fprintf(log,"EEDF-data will be recalculated!\n");
        }
    }

	fclose(tbl);
	fclose(log);

	return err;
}
