# include "1D_MainFun.h"
# include "Gas_calc.h"
# include "Chemistry_calc.h"
# include "1Dmesh.h"

int VMax[Nmax];
int ChNRe[NRmax][3];//массив номеров реагентов хим.реакций
int ChM[Nmax][NRmax],Chem[Nmax][NRmax+3];//массив мультипликаторов реакции,массив необходимых реакций
int ChRad[NRmax];//массив номеров радиационных реакций
double dEch[NRmax];//массив дефектов энергий реакций
double te_ch = 1.0, tg_ch = 1.0, tic_ch;
int Iech,Igch;

char Rtype[NRmax][10];//массив типов хим.реакций

double Kchi[NRmax][8];//массив скоростей,констант и коэффициентов скоростей хим.реакций
double RP[Nmax][LEN+2][NRmax],RC[Nmax][LEN+2][NRmax];//,RPp[LEN][Nmax][NRmax],RCp[LEN][Nmax][NRmax];
char RName[NRmax][100];

int chem_make_react(int Nedf,char *ChemFile)//формирование общего файла хим. процессов
{
	int j,n,nr,nRt,Rt,rev,J,Nn,Nvib,Nchem,k,Ki;
	char CHstr[20],Cmt[100],RTypCmt[10][100],Tcmt[10];
	char symb[] = "------------------------------------------------------------";
	char Creac[10][10],Micmt[10][10];
	double x1[10],CKj[10];

	FILE *react;
	//react = fopen("N2_chemistry.txt", "r");
	//react = fopen("N2_chemistry_new.txt", "r");
	//react = fopen("N2_chemistry_new-2.txt", "r");
	//react = fopen("N2_chemistry_new-test.txt", "r");
	//react = fopen("N2_chemistry_new-test.txt", "r");
	//react = fopen("N2_chemistry_old.txt", "r");
	//react = fopen("test-react.txt", "r");
	react = fopen(ChemFile, "r");

	FILE *rset;
	rset = fopen("ReactionSet.txt", "a+");

	//Printing diagnosic message************************************************
	printf("Chemistry info:\n\n");

	//считывание комментария
	fscanf(react,"%s",&Cmt);
	do
	{
		fscanf(react,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	Rt = 0;
	Nchem = Nedf;
	Nvib = 0;
	fscanf(react,"%s",&Cmt);
	while(strcmp(Cmt,"END")!=0)
	{
		strcpy(RTypCmt[Rt],Cmt);//RTypCmt[Rt] = Cmt;

		//printing name of reactions' block
		fprintf(rset,"\n%s\n",symb);
		fprintf(rset,"\n%s\n",RTypCmt[Rt]);
		fprintf(rset,"\n%s\n",symb);

		fscanf(react,"%s",&Cmt);
		nRt = 0;
		fscanf(react,"%s",&CHstr);
		do//read reactions inside the block
		{
			//read i-th reaction
			//left part of the reaction
			j = 0;
			while((strcmp(CHstr,"->")!=0)&&(strcmp(CHstr,"<->")!=0))
			{
				strcpy(Creac[j],CHstr);//Creac[j] = CHstr;
				j += 1;
				fscanf(react,"%s",&CHstr);
			}

			//reverse or not?
			if(!strcmp(CHstr,"<->"))
				rev = 1;
			else
				rev = 0;

			strcpy(Creac[j],"->");//Creac[j] = "->";
			j += 1;

			//right part of the reaction
			fscanf(react,"%s",&CHstr);
			while(strcmp(CHstr,";")!=0)
			{
				strcpy(Creac[j],CHstr);//Creac[j] = CHstr;
				j += 1;
				fscanf(react,"%s",&CHstr);
			}
			J = j;

			//coefficient
			fscanf(react,"%s",&Tcmt);
            if(!strcmp(Tcmt,"Te"))
                Ki = 3;
            if(!strcmp(Tcmt,"Tgas") || !strcmp(Tcmt,"Radiat"))
                Ki = 3;
            if(!strcmp(Tcmt,"VV"))
                Ki = 3;
            if(!strcmp(Tcmt,"VV'"))
                Ki = 3;
            if(!strcmp(Tcmt,"VT"))
                Ki = 3;

            for(k=0;k<Ki;k++)
                fscanf(react,"%lf",&CKj[k]);

			//different components
			fscanf(react,"%s",&CHstr);
			n = 0;
			Nn = 0;
			if(!strcmp(CHstr,"M:"))
			{
				fscanf(react,"%s",&CHstr);
				while(strcmp(CHstr,";")!=0)
				{
					strcpy(Micmt[n],CHstr);//Micmt[n] = CHstr;
					fscanf(react,"%lf",&x1[n]);
					n += 1;

					fscanf(react,"%s",&CHstr);
				}

				fscanf(react,"%s",&CHstr);
				Nn = n-1;
			}

			//add reactions to the ReactionSet.txt
            nr = Nchem+nRt;
            if((!strcmp(Tcmt,"VV"))||(!strcmp(Tcmt,"VV'"))||(!strcmp(Tcmt,"VT")))
            {
                fclose(rset);

                Nvib = chem_VV_VT_make(nr,Tcmt,Creac[0],Creac[2],CKj,Ki);

                nRt += Nvib;
            }
            else
            {
                for(n=0;n<=Nn;n++)
                {
                    nr += 1;
                    fprintf(rset,"R%d\t",nr);
                    for(j=0;j<J;j++)
                    {
                        if((!strcmp(Creac[j],"M"))&&(Nn!=0))
                            fprintf(rset,"%s\t",Micmt[n]);
                        else
                            fprintf(rset,"%s\t",Creac[j]);
                    }
                    fprintf(rset,";\t");

                    if(Nn==0)
                        fprintf(rset,"%s\t%.3e\t",Tcmt,CKj[0]);
                    else
                        fprintf(rset,"%s\t%.3e\t",Tcmt,CKj[0]*x1[n]);
                    for(k=1;k<Ki;k++)
                        fprintf(rset,"%.3lf\t",CKj[k]);
                    fprintf(rset,"%s\n",CHstr);
                    nRt += 1;

                    if(rev==1)
                    {
                        nr += 1;
                        fprintf(rset,"R%d\t",nr);
                        for(j=J-1;j>=0;j--)
                        {
                            if((!strcmp(Creac[j],"M"))&&(Nn!=0))
                                fprintf(rset,"%s\t",Micmt[n]);
                            else
                            fprintf(rset,"%s\t",Creac[j]);
                        }
                        fprintf(rset,";\t");
                        fprintf(rset,"Rev\t\t%s\n",CHstr);
                        nRt += 1;
                    }
                }
			}

			//!!***************************
			//введение исключения чтения
			//!!***************************

			fscanf(react,"%s",&Cmt);
			if(strcmp(Cmt,symb))
				strcpy(CHstr,Cmt);

		}while(strcmp(Cmt,symb)!=0);

		printf("%d\t%s were added to ReactionSet.txt\n",nRt,RTypCmt[Rt]);
		Rt += 1;
		Nchem += nRt;

		fscanf(react,"%s",&Cmt);
	}
	printf("Total of %d reactions were added to ReactionSet.txt\n\n%s\n\n",Nchem,symb);

	fclose(react);
	fclose(rset);

	return Nchem;
}
int chem_VV_VT_make(int Nchem,char Cmt[],char VSpec1[],char VSpec2[],double CKi[],int Ki)//формирование VV_VT_реакций
{
    int Nvib = 0,v1,v2,Vmax1,Vmax2,n,n1,n2,i,quant;
    char str[10],str1[10];

    FILE *rset;
	rset = fopen("ReactionSet.txt", "a+");

    //определение 1-го компонента и максимума колебательного номера
    for(n=0;n<Nmax;n++)
    {
        if(!strcmp(VSpec1,Spec[n]))
        {
            n1 = n;
            if(!strcmp(Cmt,"VT"))
                n1 = n-1;

            i = 0;
            memset(str, 0, sizeof(str));
            memset(str1, 0, sizeof(str1));
            while(VSpec1[i]!='(')
            {
                str[i] = VSpec1[i];
                i++;
            }

            Vmax1 = 0;
            do
            {
                Vmax1 += 1;
                sprintf(str1,"%s(%d)",str,Vmax1);
            }while(!strcmp(str1,Spec[n1+Vmax1]));

            VMax[n1] = Vmax1;

            break;
        }
    }

    //определение 2-го компонента
    for(n=0;n<Nmax;n++)
    {
        if(!strcmp(VSpec2,Spec[n]))
        {
            n2 = n;
            break;
        }
    }

	if(!strcmp(Cmt,"VV"))
    {
        quant = n2 - n1;

        for(v1=0;v1<Vmax1-quant;v1++)//(v1=0;v1<Vmax1-1;v1++)
        {
            for(v2=quant;v2<Vmax1;v2++)//(v2=1;v2<Vmax1;v2++)
            {
                fprintf(rset,"R%d\t %s\t +\t %s\t->\t %s\t +\t %s\t;\t",Nchem+1+Nvib,Spec[n1+v1],Spec[n1+v2],Spec[n1+v1+1],Spec[n1+v2-1]);
                if(Nvib==0)
                {
                    fprintf(rset,"VV\t %.2e\t",CKi[0]);
                    for(i=1;i<Ki;i++)
                        fprintf(rset,"%.2lf\t",CKi[i]);
                }
                else
                    fprintf(rset,"None\t");
                fprintf(rset,"//ref\n");

                Nvib += 1;
            }
        }
    }

    if(!strcmp(Cmt,"VV'"))
    {
        i = 0;
        memset(str, 0, sizeof(str));
        memset(str1, 0, sizeof(str1));
        while(VSpec2[i]!='(')
        {
            str[i] = VSpec2[i];
            i++;
        }

        n2 -= 1;
        Vmax2 = 0;
        do
        {
            Vmax2 += 1;
            sprintf(str1,"%s(%d)",str,Vmax2);
        }while(!strcmp(str1,Spec[n2+Vmax2]));

        VMax[n2] = Vmax2;

        for(v1=0;v1<Vmax1-1;v1++)
        {
            for(v2=1;v2<Vmax2;v2++)
            {
                fprintf(rset,"R%d\t %s\t +\t %s\t->\t %s\t +\t %s\t;\t",Nchem+1+Nvib,Spec[n1+v1],Spec[n2+v2],Spec[n1+v1+1],Spec[n2+v2-1]);
                if(Nvib==0)
                {
                    fprintf(rset,"VV'\t %.2e\t",CKi[0]);
                    for(i=1;i<Ki;i++)
                        fprintf(rset,"%.2lf\t",CKi[i]);
                }
                else
                    fprintf(rset,"None\t");
                fprintf(rset,"//ref\n");

                Nvib += 1;

                fprintf(rset,"R%d\t %s\t +\t %s\t->\t %s\t +\t %s\t;\tNone\t",Nchem+1+Nvib,Spec[n1+v1+1],Spec[n2+v2-1],Spec[n1+v1],Spec[n2+v2]);
                fprintf(rset,"//ref\n");

                Nvib += 1;
            }
        }
    }

    if(!strcmp(Cmt,"VT"))
    {
        for(v1=1;v1<Vmax1;v1++)
        {
            fprintf(rset,"R%d\t %s\t +\t %s\t->\t %s\t +\t %s\t;\t",Nchem+1+Nvib,Spec[n1+v1],Spec[n2],Spec[n1+v1-1],Spec[n2]);
            if(Nvib==0)
            {
                fprintf(rset,"VT\t %.2e\t",CKi[0]);
                for(i=1;i<Ki;i++)
                    fprintf(rset,"%.2lf\t",CKi[i]);
            }
            else
                fprintf(rset,"None\t");
            fprintf(rset,"//ref\n");

            Nvib += 1;

            fprintf(rset,"R%d\t %s\t +\t %s\t->\t %s\t +\t %s\t;\tNone\t",Nchem+1+Nvib,Spec[n1+v1-1],Spec[n2],Spec[n1+v1],Spec[n2]);
            fprintf(rset,"//ref\n");

            Nvib += 1;
        }
    }

    fclose(rset);

    return Nvib;
}
void chem_read_react(int Nchem,int N)//считывание процессов
{
	int i,II,j,k,n,err,S_err,Num,n_rad;
	char Cmt[300],Rstr[10];
	char symb[] = "------------------------------------------------------------";

	FILE *react;
	react = fopen("ReactionSet.txt", "r");

	//считывание комментария
	fscanf(react,"%s",&Cmt);
	do
	{
		fscanf(react,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	//Reading of the reactions**********************************************************************
	Num = 0;
	S_err = 0;
	n_rad = 0;
	for(i=0;i<Nchem;i++)
	{
		if(i==0)
			fscanf(react,"%s",&Cmt);//номер реакции

		memset(RName[i], 0, sizeof(RName[i]));
		sprintf(RName[i],"(R%d)",i+1);

		//считывание левой части реакции+формирование матрицы левых частей
		j = 0;
        fscanf(react,"%s",&Rstr);
		while(strcmp(Rstr,"->")!=0)
		{
			strcat(RName[i],Rstr);

			if(!strcmp(Rstr,"M"))
			{
				ChNRe[i][j] = -1;
				j += 1;
			}
			else if(!strcmp(Rstr,"+"))
			{}
			else
			{
				err = 0;
				for(n=0;n<N;n++)
				{
					if(!strcmp(Rstr,Spec[n]))
					{
						ChNRe[i][j] = n+1;
						ChM[n][i] -= 1;
						err = 1;
						break;
					}
				}
				if(err==0)
				{
					S_err += 1;
					printf("Unexpected specie: Left - R#=%d\t %s\n",i+1,Rstr);
				}

				j += 1;
			}
			fscanf(react,"%s",&Rstr);
		}

        strcat(RName[i],"->");

		//считывание правой части реакции
		fscanf(react,"%s",&Rstr);
		while(strcmp(Rstr,";")!=0)
		{
			strcat(RName[i],Rstr);

			if(!strcmp(Rstr,"M"))
			{}
			else if(!strcmp(Rstr,"+"))
			{}
			else
			{
				err = 0;
				for(n=0;n<N;n++)
				{
					if(!strcmp(Rstr,Spec[n]))
					{
						ChM[n][i] += 1;
						err = 1;
						break;
					}
				}
				if(err==0)
				{
					S_err += 1;
					printf("Unexpected specie: Right - R#=%d\t %s\n",i+1,Rstr);
				}
			}

			fscanf(react,"%s",&Rstr);
		}

		fscanf(react,"%s",&Rtype[i]);//тип реакции
		if(!strcmp(Rtype[i],"EEDF"))
			II = 0;
		if(!strcmp(Rtype[i],"Te"))
			II = 3;
		if(!strcmp(Rtype[i],"Tgas"))
			II = 3;
        if(!strcmp(Rtype[i],"Radiat"))
        {
            II = 3;

            n_rad += 1;
            ChRad[0] = n_rad;
            ChRad[n_rad] = i;
        }
		if(!strcmp(Rtype[i],"Tgas-VT"))
			II = 3;
        if((!strcmp(Rtype[i],"VV"))||(!strcmp(Rtype[i],"VV'"))||(!strcmp(Rtype[i],"VT")))
            II = 3;
		if((!strcmp(Rtype[i],"Rev"))||(!strcmp(Rtype[i],"None")))
			II = 0;

		//считывание коэф-та скорости реакции
		for(k=0;k<II;k++)
			fscanf(react,"%lf",&Kchi[i][k]);
		fscanf(react,"%s",&Cmt);

		Num += 1;

		//считывание комментария
		fscanf(react,"%s",&Cmt);

		if(!strcmp(Cmt,symb))
		{
			do
			{
				fscanf(react,"%s",&Cmt);
			}while(strcmp(Cmt,symb)!=0);

			fscanf(react,"%s",&Cmt);
		}

	}
	fclose(react);

	//печать сообщения об ошибках или успехе
	if((S_err==0) && (Num==Nchem))
		printf("All reactions were carefully read\n\n");
	else
	{
		char a;

		printf("!!! Please, upgrade initial file!!Mistakes Sum = %d!!!\n",S_err);
		//printf("!!! Reactions read = %d, but written = %d!!!\n",Num,Nchem);

        scanf("%c",&a);
        exit(1);
	}

	printf("%s\n",symb);

	///Расчёт_дефектов_энергии_реакции*****************************
	gas_HCpSi_calc(300,N);
	for(i=0;i<Nchem;i++)
    {
        dEch[i] = 0.0;
        for(n=1;n<N;n++)
        {
            if(ChM[n][i]!=0)
                dEch[i] += -ChM[n][i]*HCpSi[0][n]*Mi[n];///[erg]//CXi[n][0][7]*1.6e-12;
        }
    }

	///Selection of necessary reactions****************************************************
	int I,Jr,Jc,nR,m;
	FILE *rrp;
	FILE *rcp;
	char str1[20];

	for(n=0;n<N;n++)
	{
		I = 0;
		Jr = 0;
		Jc = 0;
		for(i=0;i<Nchem;i++)
		{
			if(ChM[n][i]!=0)
			{
				Chem[n][I] = i;
				I += 1;

				if(ChM[n][i]>0)
					Jr += 1;
				else
					Jc += 1;
			}
		}
		Chem[n][Nchem] = I;//число реакций, в которых участвует компонент
		Chem[n][Nchem+1] = Jr;//число реакций(приход), в которых участвует компонент
		Chem[n][Nchem+2] = Jc;//число реакций(расход), в которых участвует компонент


		//Создание файлов скоростей интересующих компонент смеси**************************
		for(nR=0;nR<NR;nR++)
		{
			if(!strcmp(Spec[n],Spec_R[nR]))
			{
				//sprintf(str1, "%s%s", &Spec_R[nR], "_R+.txt");
				//rr = fopen(str1, "a+");
				sprintf(str1, "%s%s", Spec_R[nR], "_R(t)+.txt");
				rrp = fopen(str1, "w");
				//sprintf(str1, "%s%s", &Spec_R[nR], "_R-.txt");
				//rc = fopen(str1, "a+");
				sprintf(str1, "%s%s", Spec_R[nR], "_R(t)-.txt");
				rcp = fopen(str1, "w");

				//fprintf(rr,"t,s\t");
				fprintf(rrp,"tch,s\t");
				//fprintf(rc,"t,s\t");
				fprintf(rcp,"tch,s\t");

				for(j=0;j<I;j++)
				{
					m = Chem[n][j];

					if(ChM[n][m]<0)
					{
					//	fprintf(rc,"R%d\t",m+1);
						//fprintf(rcp,"R%d\t",m+1);
						fprintf(rcp,"%s\t",RName[m]);
					}
					else
					{
					//	fprintf(rr,"R%d\t",m+1);
						//fprintf(rrp,"R%d\t",m+1);
						fprintf(rrp,"%s\t",RName[m]);
					}
				}
				//fprintf(rr,"SUM\n");
				fprintf(rrp,"RSUM\n");
				//fprintf(rc,"SUM\n");
				fprintf(rcp,"RSUM\n");

				//fclose(rr);
				fclose(rrp);
				//fclose(rc);
				fclose(rcp);

				break;
			}
		}

	}

	//Logging*******************************************************************************
	FILE *log;
	log = fopen("Log_Chemistry.txt", "w");

	//запись компонент
	fprintf(log,"Components of the mixture:\n");
	for(n=0;n<N;n++)
		fprintf(log,"%d\t%s\n",n+1,Spec[n]);
	fprintf(log,"\n");

	//запись матрицы мультипликаторов реакции
	fprintf(log,"Matrix of multipliers:\n");
	for(n=0;n<N;n++)
	{
		fprintf(log,"%s\t\t",Spec[n]);
		for(i=0;i<Nchem;i++)
			fprintf(log,"%d\t",ChM[n][i]);
		fprintf(log,"\n");
	}
	fprintf(log,"\n");

	//запись матрицы необходимых реакций
	fprintf(log,"Matrix of necessary reactions:\n");
	for(n=0;n<N;n++)
	{
		fprintf(log,"%5s\t",Spec[n]);
		for(i=0;i<=Nchem+2;i++)
			fprintf(log,"%d\t",Chem[n][i]);
		fprintf(log,"\n");
	}
	fprintf(log,"\n");

	fclose(log);
}
void chem_const(double *Kch,int Nedf,int Nchem,int N,double Tel,double Tch,double tic)// расчёт констант скоростей реакций
{
	int j,n,Sum;
	double Kf;
	double Nv,kv10;

	int xRev=0;//счётчик запуска счета констант реакций с ID - Rev
    Tel *= eV_K;//[K]

	for(j=Nedf;j<Nchem;j++)
	{
		if(!strcmp(Rtype[j],"Te"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T), k-[см^3/с], Ea-[K]
		{
			Kf = Kchi[j][0]*pow(Tel,Kchi[j][1])*exp(-Kchi[j][2]/Tel);

			Kch[j]=Kf;
		}
		if(!strcmp(Rtype[j],"Tgas") || !strcmp(Rtype[j],"Radiat"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T), k-[см^3/с], Ea-[K]
		{
			Kf = Kchi[j][0]*pow(Tch,Kchi[j][1])*exp(-Kchi[j][2]/Tch);

			Kch[j]=Kf;
		}
		if(!strcmp(Rtype[j],"Tgas-VT"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T^0.33), k-[см^3/с], Ea-[K]
		{
			double C1;

			C1 = pow(Tch,0.33);

			Kf = Kchi[j][0]*pow(Tch,Kchi[j][1])*exp(-Kchi[j][2]/C1);

			Kch[j]=Kf;
		}
        if(!strcmp(Rtype[j],"VV"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T^0.33), k-[см^3/с], Ea-[K]
        {
            kv10 = Kchi[j][0]*pow(Tch/300.0,Kchi[j][1])*exp(-Kchi[j][2]/Tch);

            Nv = chem_VV_VT_const(Kch,N,j,"VV",kv10,Tch);

            j += Nv-1;
        }
        if(!strcmp(Rtype[j],"VV'"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T^0.33), k-[см^3/с], Ea-[K]
        {
            kv10 = Kchi[j][0]*pow(Tch/300.0,Kchi[j][1])*exp(-Kchi[j][2]/Tch);

            Nv = chem_VV_VT_const(Kch,N,j,"VV'",kv10,Tch);

            j += Nv-1;
        }
        if(!strcmp(Rtype[j],"VT"))//запись констант реакций в виде - k*(T)^n*exp(-Ea/T^0.33), k-[см^3/с], Ea-[K]
        {
            kv10 = Kchi[j][0]*pow(Tch,Kchi[j][1])*exp(-Kchi[j][2]/Tch);

            Nv = chem_VV_VT_const(Kch,N,j,"VT",kv10,Tch);

            j += Nv-1;
        }
		if(!strcmp(Rtype[j],"Rev"))//вычисление констант обратных реакций
		{
			double dH,dS,Kb,Keq,ChMm;
			double Patm = 760*p0;

			if(xRev==0)
				gas_HCpSi_calc(Tch,N);

			xRev += 1;

			dS = 0;
			dH = 0;
			for(n=0;n<N;n++)
			{
				ChMm = ChM[n][j-1]*(Mi[n]*ma);///!!!!!Mass!!!

				dH += HCpSi[0][n]*ChMm;//[эрг]*1e-10*Na;/[кДж/моль]-для проверки
				dS += HCpSi[2][n]*ChMm;//[эрг/К]*1e-7*Na;//[Дж/К*моль]-для проверки
			}

			Sum = 0;
			for(n=0;n<N;n++)
				Sum += ChM[n][j-1];

			Keq = exp(dS/kb-dH/(kb*Tch))*pow(Patm/(kb*Tch),Sum);

			Kb = Kf/Keq;

			Kch[j] = Kb;
		}

        //printf("R#%d\tKch=%.2e[cm^-3]\n",j+1,Kch[j]);

	}

	//Logging_Kch***********************************************************
	if(tic==0.0)//(dot==Ndots)
	{
		FILE *log;
		if(tic==0)
            log = fopen("Log_Kch.txt", "w");
        else
            log = fopen("Log_Kch.txt", "a+");

        fprintf(log,"(R#)\tReaction\tt=%.2e[s]\tdE,eV\n",tic);
        for(j=0; j<Nchem; j++)
            fprintf(log,"%s\t%.2e\t%.2lf\n",RName[j],Kch[j],dEch[j]/1.6e-12);
        fprintf(log,"\n\n");

		/*else
        {
            log = fopen("Log_Kch.txt", "r+");
            //rewind(log);

            char cmt[10];
            while(strcmp(cmt,";")!=0)
                fscanf(log,"%s",&cmt);
            fprintf(log,"\tt=%.2e [s] ;",dt*nt);

            cmt[0] = '\0';
            for(j=0;j<Nchem;j++)
            {
                while(strcmp(cmt,";")!=0)
                    fscanf(log,"%s",&cmt);

                fprintf(log,"\t%.2e ;",Kch[j]);
                cmt[0] = '\0';
            }
            fprintf(log,"\n\n");
        }*/

		fclose(log);
	}

}
int chem_VV_VT_const(double *Kch,int N,int j,char vin[],double kv10,double Tch)// расчёт констант скоростей колебательных реакций
{
    double dvv,Kvv,dvt,Kvt;
    double dEv;
    int Nvib = 0;
    int n,n1,n2,quant,Vmax,Vmax1,Vmax2,v1,v2;

    gas_HCpSi_calc(Tch,N);

    n1 = ChNRe[j][0]-1;
    n2 = ChNRe[j][1]-1;

    //FILE *log;
    //log = fopen("Log_VV_VT.txt", "r+");

    if(!strcmp(vin,"VV"))
    {
        n = n1;
        quant = n2-n1;

        //Gordiets et. al_"Kineticheskie proc v gasah"_1980_Nauka
        //Billing
        //Валянский_КвантЭл_1984

        dvv=-6.85/sqrt(Tch);

        Vmax = VMax[n];
        for(v1=0;v1<Vmax-quant;v1++)//(v1=0;v1<Vmax-1;v1++)
        {
            for(v2=quant;v2<Vmax;v2++)//(v2=1;v2<Vmax;v2++)
            {
                dEv = ((HCpSi[0][n+v1]-HCpSi[0][n+v1+quant])+(HCpSi[0][n+v2]-HCpSi[0][n+v2-quant]))*Mi[n]/kb;//energy_defect[K]

                Kvv = kv10*(v1+quant)*v2*exp(dvv*abs(v1-(v2-quant)))*(1.5-0.5*exp(dvv*abs(v1-(v2-quant))));//*exp(dEv*(v1-(v2-quant))/Tch);

                if(dEv<0)//need to check!!!!!!
                    Kvv *= exp(dEv/Tch);

                Kch[j+Nvib] = Kvv;

                Nvib++;

                //fprintf(log,"R#%d\t%.2lf\n",j+Nvib-1,dEv);
            }
        }

    }
    if(!strcmp(vin,"VV'"))
    {
        n2 -= 1;//
        Vmax1 = VMax[n1];
        Vmax2 = VMax[n2];

        //in analogy with VV: MUST be validated!!!
        //Gordiets et. al_"Kineticheskie proc v gasah"_1980_Nauka
        dvv=-6.85/sqrt(Tch);

        for(v1=0;v1<Vmax1-1;v1++)
        {
            for(v2=1;v2<Vmax2;v2++)
            {
                //energy_defect[K]
                dEv = (HCpSi[0][n1+v1]*Mi[n1+v1]-HCpSi[0][n1+v1+1]*Mi[n1+v1+1])/kb;
                dEv += (HCpSi[0][n2+v2]*Mi[n2+v2]-HCpSi[0][n2+v2-1]*Mi[n2+v2-1])/kb;

                Kvv = kv10*(v1+1)*v2*exp(dvv*abs(v1-(v2-1)))*(1.5-0.5*exp(dvv*abs(v1-(v2-1))));//*exp(dEv*(v1-(v2-1))/Tch);

                if(dEv<0)//need to check!!!!!!
                    Kvv = Kvv*exp(dEv/Tch);

                Kch[j+Nvib] = Kvv;
                Kch[j+Nvib+1] = Kvv*exp(-dEv/Tch);

                Nvib += 2;
            }
        }
    }
    if(!strcmp(vin,"VT"))
    {
        n1 -= 1;
        Vmax = VMax[n1];

        //Approximation!!!

        if(!strcmp(Spec[n1],"N2(0)"))
        {

            if(!strcmp(Spec[n2],"N2(0)"))
            {
                dvt = 2.87/pow(Tch,0.333333);//Billing

                for(v1=1;v1<Vmax;v1++)
                {
                    Kvt = kv10*v1*exp((v1-1)*dvt);
                    dEv = (HCpSi[0][n1+v1]-HCpSi[0][n1+v1-1])*Mi[n1]/kb;//[K]

                    Kch[j+Nvib] = Kvt;
                    Kch[j+Nvib+1] = Kvt*exp(-dEv/Tch);

                    //fprintf(log,"R#%d\t%.2lf\n",j+Nvib,dEv);
                    //fprintf(log,"R#%d\t%.2lf\n",j+Nvib+1,-dEv);

                    Nvib += 2;
                }

            }
            if(!strcmp(Spec[n2],"N"))//Capitelli,Ferreira_"Plasma_Kinetics"_2000
            {
                double Ev, beta = 0.065;//Capitelli_2000 //beta = 0.065;//Capitelli_2006
                //int dvmax = 5;

                for(v1=1;v1<Vmax;v1++)
                {
                    Ev = (HCpSi[0][n1+v1]-HCpSi[0][n1])*Mi[n1]/kb;//[K]
                    Kvt = kv10*exp(beta*Ev/Tch);

                    /*
                    //for multi-quantum transitions
                    if(v1<=dvmax)
                        Kvt = Kvt/v1;
                    else
                        Kvt = Kvt/dvmax;
                    */

                    dEv = (HCpSi[0][n1+v1]-HCpSi[0][n1+v1-1])*Mi[n1]/kb;//[K]

                    Kch[j+Nvib] = Kvt;
                    Kch[j+Nvib+1] = Kvt*exp(-dEv/Tch);

                    Nvib += 2;
                }
            }
        }
    }

    return Nvib;
}
void chem_runge_kutta4(int I,double *Ni,int Nion,int N,double *Kch,double *Sch,int Nchem,double *Wrad,double dt,double dte,double *tg_chem,double *te_chem,double tic,int dot)//расчёт по времени методом Рунге-Кутта 4 порядка
{
	int i,n,m,s,j,l;
	double al[4] = {0,dt/2,dt/2,dt};
	double Ci[N];
	double Rch[Nchem],Rch_rk[N][4],RRg;
	double M,Res,Yi,dEi;

    for(n=0;n<N;n++)
        Ci[n] = *(Ni+n*(LEN+2));

	///**********************расчёт по явной схеме(Рунге-Кутта_4)**************************
	for(i=0;i<1;i++)
	//for(i=0;i<4;i++)
	{
		if(i>0)
        {
            for(n=0;n<N;n++)
                Ci[n] = *(Ni+n*(LEN+2))+al[i]*Rch_rk[n][i-1];
        }

		///Расчет скоростей химических реакций************************************************************
		M = 0.0;
		for(n=1;n<N;n++)
			M += Ci[n];

		for(j=0;j<Nchem;j++)
		{
			Res = 1.0;
			for(l=0;l<3;l++)
			{
				if(ChNRe[j][l]==-1)
					Yi = M;
				else if(ChNRe[j][l]==0)
					Yi = 1.0;
				else
				{
					n = ChNRe[j][l]-1;
					Yi = Ci[n];
				}
				Res *= Yi;
			}

			Rch[j] = Kch[j]*Res;
		}

        ///Расчет Радиационных потерь тепла******************************************************************
        if(i==0)
        {
            *Wrad = 0.0;
            for(s=1;s<ChRad[0];s++)
            {
                m = ChRad[s];
                *Wrad += Rch[m]*dEch[m];///[erg/cm3/s]
            }
        }

		///Коэф-ты к правой части реакций(necessary processes)************************************************
		int Num,nnR,nR=0;
        int Jr,Jc,sNR;
        double SR_prod,SR_cons;
		for(n=0;n<N;n++)
		{
		    sNR = 0;
            if((i==0) && (tic/tic_ch)>3.0  && (nnR<NR))//&& (I==1)
            {
                for(nR=0;nR<NR;nR++)
                {
                    if(!strcmp(Spec[n],Spec_R[nR]))
                    {
                        Jr=0;Jc=0;SR_prod=0;SR_cons=0;
                        sNR=1;
                        break;
                    }
                }
            }

			Num = Chem[n][Nchem];
			Rch_rk[n][i] = 0.0;
			for(j=0;j<Num;j++)
			{
				m = Chem[n][j];
				RRg = ChM[n][m]*Rch[m];
				Rch_rk[n][i] += RRg;

                ///Поиск_минимального_химического_шага*******************************************************
                if((i==0) && (RRg!=0.0) && (Ci[n]>0.0))
                {
                    if(n<=Nion)
                        te_ch = fmin(fabs(Ci[n]/RRg),te_ch);///for_charged
                    else //if(n>Nion)
                        tg_ch = fmin(fabs(Ci[n]/RRg),tg_ch);///for_neutrals
                }

                ///Вычисление долей и скоростей реакций для компонент газа***************************************
                if(sNR!=0)
                {
                    if(ChM[n][m]<0)
                    {
                        RC[nR][I][Jc] = fabs(RRg);

                        Jc += 1;
                        SR_cons += fabs(RRg);
                    }
                    else
                    {
                        RP[nR][I][Jr] = RRg;

                        Jr += 1;
                        SR_prod += RRg;
                    }

                    if(j==Num-1)
                    {
                        RP[nR][I][Jr] = SR_prod;
                        RC[nR][I][Jc] = SR_cons;

                        nnR++;
                    }
                }
            }
		}
	}

	Igch += 1;
	if(Igch==LEN)
    {
        *tg_chem = tg_ch;
        tg_ch = 1.0;
        Igch = 0;
    }

    Iech += 1;
	if(Iech==LEN)
    {
        *te_chem = te_ch;
        te_ch = 1.0;
        Iech = 0;
    }

	///обновление концентраций**********************************************************
	for(n=0;n<N;n++)
	{
        //*(Sch+n*(LEN+2) = (Rch_rk[n][0]+2*Rch_rk[n][1]+2*Rch_rk[n][2]+Rch_rk[n][3])/6;
        *(Sch+n*(LEN+2)) = Rch_rk[n][0];
    }

	///Writing_Chem-Contributions******************************************************
	if(I==LEN && (tic/tic_ch)>3.0)
    {
        chem_spec_contrib(N,Nchem,1.0,tic);
        tic_ch = tic;
    }
	///******************************************************************************

}
void chem_half_implicit(double *Ni,int N0,int N,double *Kch,int Nchem,double *Wrad,double dt,double *tg_chem,double tic,int dot)//расчёт по времени неявным методом
{
	double Y0[N],Y[N];
	double Si[N],Li[N],Ki[N],Res[N],yRes;
	double Q,A,B,C,M,Y_M,Rch;

	int n,m,j,l,L,s,NJ;

	M = 0.0;
	for(n=0;n<N;n++)
    {
        Y0[n] = *(Ni+n*(LEN+2));
        M += Y0[n];
    }

	/*
	ChNRe[j][l] - //матрица реагентов
	Chem[m][j] - //матрица номеров реакций с компонентом m
	ChM[m][n] - //матрица итога реакции с компонентом m
	*/

	for(n=N0;n<N;n++)
	{
		Li[n] = 0.0;
		Ki[n] = 0.0;
		Si[n] = 0.0;
	}

	for(j=0;j<Nchem;j++)
	{
		Y_M = 1.0;
		L = 0;
		for(l=0;l<3;l++)
		{
			if(ChNRe[j][l]==-1)
				Y_M *= M;
			if(ChNRe[j][l]>0)
			{
				m = ChNRe[j][l]-1;

				for(n=N0;n<N;n++)
				{
					if(L==0)
						Res[n] = 1.0;

					if(m!=n)
						Res[n] *= Y0[m];
				}

				L += 1;
			}
		}

		for(n=N0;n<N;n++)
		{
			Res[n] *= Kch[j]*Y_M;

			NJ = ChM[n][j];

			if(NJ==-1)
				Li[n] += Res[n];
			else if(NJ==-2)
				Ki[n] += 2.0*Res[n];
			else if(NJ>0)
				Si[n] += NJ*Res[n];
		}
	}

    for(n=N0;n<N;n++)
    {
        //Q = dt*fabs(Si[n]-(Li[n]+Ki[n]*Y0[n])*Y0[n]);
		Q = Si[n]-(Li[n]+Ki[n]*Y0[n])*Y0[n];

		tg_ch = fmin(Y0[n]/Si[n],1.0/Li[n]);
		if(Ki[n]*Y0[n]!=0)
            tg_ch = fmin(tg_ch,1.0/(Ki[n]*Y0[n]));

        //if(fabs(Y0[n]/Q)>0.5*dt)
        if(tg_ch<dt || fabs(Q)*dt>=0.05*Y0[n])
        {
            ///evential
            yRes = (Y0[n]+dt*Si[n])/(1.0+(Li[n]+Ki[n]*Y0[n])*dt);
            if(yRes<0.0 || fabs(1-yRes/Y0[n])>0.05)
            {
                yRes = chem_sor_iterator(5,1.25,Y0[n],Si[n],Li[n],Ki[n],dt,1.0e-3);
                if(yRes<0.0)
                {
                    yRes = chem_newton_iterator(10,yRes,Si[n],Li[n],Ki[n],dt,1.0e-3);
                    if(yRes<0.0 && Q>0.0)
                        yRes = Y0[n]*1.01;
                    else if(yRes<0.0 && Q<0.0)
                        yRes = Y0[n]*0.99;
                }

            }

            Y[n] = yRes;
        }
        else
			Y[n] = Y0[n]+dt*Q;

        if(Y[n]<1.e-30)
			Y[n] = 0.0;

        *(Ni+n*(LEN+2)) = Y[n];
    }

    ///Поиск_минимального_химического_шага*******************************************************
    /*if((i==0) && (RRg!=0.0) && (Ci[n]>0.0))// && (n==51))
        tg_ch = fmin(fabs(Ci[n]/RRg),tg_ch);

	Igch += 1;
	if(Igch==LEN)
    {
        *tg_chem = tg_ch;
        tg_ch = 1.0;
        Igch = 0;
    }*/

    *tg_chem = 1.0e-8;

    ///Расчет Радиационных потерь тепла**********************************************************

    *Wrad = 0.0;
    for(s=1;s<ChRad[0];s++)
    {
        j = ChRad[s];

        Rch = Kch[j];
        for(l=0;l<3;l++)
		{
			if(ChNRe[j][l]==-1)
				Y_M = M;
			else if(ChNRe[j][l]==0)
				Y_M = 1.0;
			else
			{
				n = ChNRe[j][l]-1;
				Y_M = Y[n];//Y0[n];
			}
			Rch *= Y_M;
		}

        *Wrad += Rch*dEch[m];///[erg/cm3/s]
    }
}
void chem_half_modimplicit(double *Ni,int N0,int N,double *Kch,int Nchem,double *Wrad,double dt,double *tg_chem,double tic,int dot)//расчёт по времени неявным методом
{
	double Y0[N],Y[N];
	double Si[N],Li[N],Ki[N],Res[N],yRes;
	double Q,A,B,C,M,Y_M,Rch;
	int unc,nUnc[N+1];

	int i,n,m,j,l,L,s,NJ;

	M = 0.0;
	nUnc[N] = 0;
	for(n=0;n<N;n++)
    {
        Y0[n] = *(Ni+n*(LEN+2));
        M += Y0[n];
    }

    for(n=N0;n<N;n++)
    {
        nUnc[n] = n;
        nUnc[N]++;
    }
	/*
	ChNRe[j][l] - //матрица реагентов
	Chem[m][j] - //матрица номеров реакций с компонентом m
	ChM[m][n] - //матрица итога реакции с компонентом m
	*/

	s = 0;
	do
    {
        //for(n=N0;n<N;n++)
        for(i=0;i<nUnc[N];i++)
        {
            n=nUnc[i];

            Li[n] = 0.0;
            Ki[n] = 0.0;
            Si[n] = 0.0;
        }

        for(j=0;j<Nchem;j++)
        {
            Y_M = 1.0;
            L = 0;
            for(l=0;l<3;l++)
            {
                if(ChNRe[j][l]==-1)
                    Y_M *= M;
                if(ChNRe[j][l]>0)
                {
                    m = ChNRe[j][l]-1;

                    for(n=N0;n<N;n++)
                    {
                        if(L==0)
                            Res[n] = 1.0;

                        if(m!=n)
                            Res[n] *= Y0[m];
                    }

                    L += 1;
                }
            }

            //for(n=N0;n<N;n++)
            for(i=0;i<nUnc[N];i++)
            {
                n=nUnc[i];

                NJ = ChM[n][j];
                if(NJ!=0)
                {
                    Res[n] *= Kch[j]*Y_M;

                    if(NJ==-1)
                        Li[n] += Res[n];
                    else if(NJ==-2)
                        Ki[n] += 2.0*Res[n];
                    else if(NJ>0)
                        Si[n] += NJ*Res[n];
                }
            }
        }

        //for(n=N0;n<N;n++)
        for(i=0;i<nUnc[N];i++)
        {
            n=nUnc[i];

            Q = Si[n]-(Li[n]+Ki[n]*Y0[n])*Y0[n];

            tg_ch = fmin(Y0[n]/Si[n],1.0/Li[n]);
            if(Ki[n]*Y0[n]!=0)
                tg_ch = fmin(tg_ch,1.0/(Ki[n]*Y0[n]));

            if(tg_ch<dt)
            {
                if(s<5)
                    yRes = chem_sor_iterator(1,1.25,Y0[n],Si[n],Li[n],Ki[n],dt,1.0e-3);
                else
                    yRes = chem_newton_iterator(1,Y0[n],Si[n],Li[n],Ki[n],dt,1.0e-3);

                Y[n] = yRes;
            }
            else
                Y[n] = Y0[n]+dt*Q;

            if(Y[n]<1.e-30)
                Y[n] = 0.0;

            if(fabs(1-Y[n]/Y0[n])>1.0e-3)///no_converge
            {
                nUnc[unc] = n;
                unc++;
            }

             Y0[n] = Y[n];
        }
        nUnc[N] = unc;
        unc = 0;

        s++;

    }while(nUnc[N]!=0 || s<25);

    for(n=N0;n<N;n++)
        *(Ni+n*(LEN+2)) = Y[n];

    /*///Поиск_минимального_химического_шага*******************************************************
    if((i==0) && (RRg!=0.0) && (Ci[n]>0.0) && (n==51))
        tg_ch = fmin(fabs(Ci[n]/RRg),tg_ch);

	Igch += 1;
	if(Igch==LEN)
    {
        *tg_chem = tg_ch;
        tg_ch = 1.0;
        Igch = 0;
    }
    */

    *tg_chem = 5.0e-7;

    ///Расчет Радиационных потерь тепла**********************************************************

    *Wrad = 0.0;
    for(s=1;s<ChRad[0];s++)
    {
        j = ChRad[s];

        Rch = Kch[j];
        for(l=0;l<3;l++)
		{
			if(ChNRe[j][l]==-1)
				Y_M = M;
			else if(ChNRe[j][l]==0)
				Y_M = 1.0;
			else
			{
				n = ChNRe[j][l]-1;
				Y_M = Y[n];//Y0[n];
			}
			Rch *= Y_M;
		}

        *Wrad += Rch*dEch[m];///[erg/cm3/s]
    }
}
void chem_spec_contrib(int N,int Nchem,double Rth,double tic)//вывод скоростей реакций для выбранных компонент
{
	FILE *rrp;
	FILE *rcp;

	char str1[20];

	int nR,i,j,n,m,k,s,Jr,Jc,Num;
	double Rr[LEN+2],Rav;
	double Rp_av[Nchem],Rc_av[Nchem];

	for(nR=0;nR<NR;nR++)
	{
	    ///Timeline_of_Rates***********************************************************************
		sprintf(str1, "%s%s", Spec_R[nR], "_R(t)+.txt");
		rrp = fopen(str1, "a+");
		sprintf(str1, "%s%s", Spec_R[nR], "_R(t)-.txt");
		rcp = fopen(str1, "a+");

		for(n=0;n<N;n++)
		{
			if(!strcmp(Spec[n],Spec_R[nR]))
			{
				Jr = Chem[n][Nchem+1];//число реакций(приход), в которых участвует компонент
				Jc = Chem[n][Nchem+2];//число реакций(расход), в которых участвует компонент

				break;
			}
		}

        fprintf(rrp,"%.2e\t",tic);
		for(j=0;j<=Jr;j++)
			fprintf(rrp,"%.2e\t",RP[nR][1][j]);
		fprintf(rrp,"\n");

		fprintf(rcp,"%.2e\t",tic);
		for(j=0;j<=Jc;j++)
			fprintf(rcp,"%.2e\t",RC[nR][1][j]);
		fprintf(rcp,"\n");

		fclose(rrp);
		fclose(rcp);

//		double Rp_av[Jr+1],Rc_av[Jc+1];

        ///Rate_Fraction_over_length****************************************************************
        for(j=0;j<=Jr;j++)///Production
        {
            for(i=1;i<=LEN;i++)
                Rr[i] = RP[nR][i][j];
            Rav = LENaverage(Rr,Geom);

            for(i=1;i<=LEN;i++)
                Rr[i] = RP[nR][i][Jr];
            Rav *= 100/LENaverage(Rr,Geom);
            Rp_av[j] = Rav;///in_percent_%
        }

        for(j=0;j<=Jc;j++)///Consumption
        {
            for(i=1;i<=LEN;i++)
                Rr[i] = RC[nR][i][j];
            Rav = LENaverage(Rr,Geom);

            for(i=1;i<=LEN;i++)
                Rr[i] = RC[nR][i][Jc];
            Rav *= 100/LENaverage(Rr,Geom);
            Rc_av[j] = Rav;///in_percent_%
        }
        ///*****************************************************************************************

		///Spatial_Rates_Distribution***************************************************************
        sprintf(str1, "%s%s", Spec_R[nR], "_RLen+.txt");
        rrp = fopen(str1, "w");

        sprintf(str1, "%s%s", Spec_R[nR], "_RLen-.txt");
        rcp = fopen(str1, "w");

        fprintf(rrp,"r_cm\t");
        fprintf(rcp,"r_cm\t");

        Num = Chem[n][Nchem];
        k=0;
        s=0;
        for(j=0;j<Num;j++)
        {
            m = Chem[n][j];

            if(ChM[n][m]>0)
            {
                if(Rp_av[k]>=Rth)
                    fprintf(rrp,"%s_\t",RName[m]);
                k++;
            }
            else
            {
                if(Rc_av[s]>=Rth)
                    fprintf(rcp,"%s_\t",RName[m]);
                s++;
            }
        }
        fprintf(rrp,"%s_RSUM+_\n",Spec_R[nR]);
        fprintf(rcp,"%s_RSUM-_\n",Spec_R[nR]);

        fprintf(rrp,"<R,%%>\t");
        for(j=0;j<=Jr;j++)
        {
            if(Rp_av[j]>=Rth)
                fprintf(rrp,"%.2lf\t",Rp_av[j]);
        }
        fprintf(rrp,"\n");

        fprintf(rcp,"<R,%%>\t");
        for(j=0;j<=Jc;j++)
        {
            if(Rc_av[j]>=Rth)
                fprintf(rcp,"%.2lf\t",Rc_av[j]);
        }
        fprintf(rcp,"\n");

        for(i=1;i<=LEN;i++)
        {
            fprintf(rrp,"%.3lf\t",0.5*(l[i+1]+l[i]));
            for(j=0;j<=Jr;j++)
            {
                if(Rp_av[j]>=Rth)
                    fprintf(rrp,"%.2e\t",RP[nR][i][j]);
            }
            fprintf(rrp,"\n");

            fprintf(rcp,"%.3lf\t",0.5*(l[i+1]+l[i]));
            for(j=0;j<=Jc;j++)
            {
                if(Rc_av[j]>=Rth)
                    fprintf(rcp,"%.2e\t",RC[nR][i][j]);
            }
            fprintf(rcp,"\n");
        }

        fclose(rrp);
        fclose(rcp);
	}
}
double chem_sor_iterator(int I,double w,double y0,double S,double L,double K,double dt,double eps)
{
    double y,A,B,Norm;
    int i;

    i = 0;
    do
    {
        A = 1.0+dt*(L+K*y0);///impl-part
        B = y0+dt*S;///expl-part

        y = y0*(1-w) + w*B/A;

        Norm = fabs(1-y/y0);

        y0 = y;

        i++;

    }while(Norm>eps || i<I);

    if(Norm>eps)
        return -1.0;
    else
        return y;
}
double chem_newton_iterator(int I,double y0,double S,double L,double K,double dt,double eps)
{
    double y,A,B,Norm,f,F;
    int i;

    i = 0;
    do
    {
        A = 1.0+dt*(L+K*y0);///impl-part
        B = y0+dt*S;///expl-part

        f = (2.0*K*y0+L)*dt;

        F = A*y0-B;

        y = y0 - F/f;

        Norm = fabs(1-y/y0);

        y0 = y;

        i++;

    }while(Norm>eps || i<I);

    if(Norm>eps)
        return -1.0;
    else
        return y;
}
