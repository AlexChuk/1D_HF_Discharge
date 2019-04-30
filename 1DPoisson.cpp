# include "1D_MainFun.h"
# include "1DPoisson.h"
# include "1Dmesh.h"
# include "Initials.h"

extern double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];
extern int vL,vR;
extern double VolL,VolR,WextL,WextR,EpsL,EpsR,DepL,DepR;

double Poisson_SORsolve(double *Fi,double *E,double Vext,double *Xi,double *Ni,double *Ngas,int Npos,int Nneg,double tic)
{
	/*
	**************************************************************

	Poisson equation solution with
	Successive over Relaxation (SOR) Method

	if w = 1.0
	SOR = Gauss-Seidel

	**************************************************************
	*/

    /*
	!!!!!!!!!Необходимые доработки:
        - Поток зарядов на стенку
        - улучшить алгоритм сходимости (не гонять весь диапазон по координате)_ для "частых" сеток
        - определить алгоритм выбора константы w
	*/

	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    int i,n;
	double Res,Ch,RHS[LEN+2];
    int cnt,Conv;//conv[LEN+2];
    double VL,VR,SigL,SigR,Norm;
    double E0[LEN+1],Fi0[LEN+2],dEf,Emax;

    ///Initial_field*********************************
    for(i=0;i<=LEN;i++)
    {
        Fi0[i] = Fi[i];
        E0[i] = E[i];
    }
    Fi0[LEN+1] = Fi[LEN+1];

    ///Volume_Charge_calculation*********************
    for(i=0;i<=LEN+1;i++)
    {
        Ch = -Ni[i] ;//electrons//Ni[i];
        for(n=1;n<=Npos;n++)//positive ions
            Ch += Ni[n*(LEN+2)+i];//Ni[n*(LEN+2)+i];

        for(n>Npos;n<=Npos+Nneg;n++)//negative ions
            Ch += -Ni[n*(LEN+2)+i] ;

        //if(fabs(Ch/Xi[i])<1.e-3)
            //Ch = 0.0;

        RHS[i] = -4*pi*e*Ch;//Ngas[i];
    }

    ///Surface_charge_density************************
    SigL = 0.0;
    SigR = 0.0;///dSig/dt = - Jrad;

    if(Vext!=0.0)///if_ExtVoltage_changes
    {
        VL = VolL;
        if(VolL!=0.0)
            VL = Vext;

        VR = VolR;
        if(VolR!=0.0)
            VR = Vext;
    }
    else///if_ExtVoltage_const
    {
        VL = VolL;
        VR = VolR;
    }
    if(WextL!=0.0)///[MHz]
        VL *= sin(2*pi*WextL*tic);
    if(WextR!=0.0)
        VR *= sin(2*pi*WextR*tic);

    ///SOR_cycle_Fi_calculation**********************
    /*cnt = 0;
    double w = 1.25;//Gauss-Seidel if w = 1.0
	do
	{
        Poisson_boundary(0,Fi,vL,VolL,EpsL,DepL,RHS[0],SigL);///Left_boundary

        Conv = 0;
        for(i=1;i<=LEN;i++)
        {
            //Res = Fi[i]*(1-w) + w*(RHS[i]-GF_L[i]*Fi[i-1]-GF_R[i]*Fi[i+1])/(-GF_C[i]);
            Res = Fi[i] + w*((RHS[i]-GF_L[i]*Fi[i-1]-GF_R[i]*Fi[i+1])/(-GF_C[i])-Fi[i]);
            ///Res = Fi[i] + w*((RHS[i]-Fi[i-1]-Fi[i+1])/(-2)-Fi[i]);///test


            //Norm = (Fi[i]-Res)*(Fi[i]-Res);
            //Norm = sqrt(Norm);
            //Norm = fabs(Res-Fi[i]);
            Norm = fabs(1-Res/Fi[i]);

            if(Norm<1.0e-3)//(fabs(Fi[i]-Res)<exact);//1.0-Res/Fi[i])<1.0e-5_i=1 &&
                Conv++;

            Fi[i] = Res;///[СГС]

            //if(Fi[i]<1.0e-5)
                //Fi[i] = 0.0;
        }

        Poisson_boundary(LEN+1,Fi,vR,VolR,EpsR,DepR,RHS[LEN+1],SigR);///Right_boundary

        cnt++;

	}while(Conv<LEN || cnt<30);//(Conv<LEN
    */

    ///SWEEP_method***************************************************************
    double  A,B,C,F,den,
            al[LEN+1],bet[LEN+1];

    //Poisson_boundary(0,Fi,vL,VL,EpsL,DepL,RHS[0],SigL);///Left_boundary
    Poisson_boundary(LEN+1,Fi,vR,VR,EpsR,DepR,RHS[LEN+1],SigR);///Right_boundary
    for(i=1;i<=LEN;i++)
    {
        ///SWEEP Coefficients:
        A = GF_L[i];///[i-1](left_cell)
        B = GF_R[i];///[i+1](right_cell)
        C = -GF_C[i];///[i](center_cell)
        F = RHS[i];///RHS-part

        if(i==1)//see SpecTransport();
        {
            //al[i-1] = 1.0;///see_SpecTransportBoundary();
            al[i-1] = 0.0;///see_SpecTransportBoundary();
            bet[i-1] = 0.0;//0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;
    }

    ///Reverse_sweep_cycle************************************************
    for(i=LEN;i>=0;i--)
        Fi[i] = al[i]*Fi[i+1]+bet[i];

    ///Potential_Smoothing???
    for(i=0;i<=LEN+1;i++)
        //Fi[i] = 0.9*Fi0[i]+0.1*Fi[i];
        Fi[i] = 0.7*Fi0[i]+0.3*Fi[i];

    ///Field_Calculation**************************************************
    for(i=0;i<=LEN;i++)
        E[i] = -2.0*(Fi[i+1]-Fi[i])/(l[i+2]-l[i]);///[СГС]

    //dEf = 0.0;
    Emax = 0.0;
    for(i=0;i<=LEN;i++)
        //dEf = fmax(fabs(E[i]-E0[i]),dEf);
        Emax = fmax(fabs(E[i]),Emax);

    return Emax;
}
void Poisson_boundary(int i,double *Fi,int v,double Vext,double Eps,double Dep,double Rhs,double Sig)
{
    int iB,iV,iV2;
    double dl,l0,l1,l2,dl2;

    if(i==0)///Left_boundary
    {
        iB = 0;///boundary
        iV = 1;///volume
        iV2 = 2;///volume2

        dl = 0.5*(l[2]-l[0]);
        l0 = 0.5*(l[1]+l[0]);
        l1 = 0.5*(l[2]+l[1]);
        l2 = 0.5*(l[3]+l[2]);
        dl2 = 0.5*(l[3]-l[1]);
    }
    else///Right_boundary
    {
        iV = LEN;///volume
        iV2 = LEN-1;///volume2
        iB = LEN+1;///volume

        dl = 0.5*(l[LEN+2]-l[LEN]);
        l0 = 0.5*(l[LEN+2]+l[LEN+1]);
        l1 = 0.5*(l[LEN+1]+l[LEN]);
        l2 = 0.5*(l[LEN]+l[LEN-1]);
        dl2 = 0.5*(l[LEN+1]-l[LEN-1]);
    }

    ///Cases:
    if(v==0)///defined
        Fi[iB] = Vext;
    if(v==1)///symmetry
        //Fi[iB] = Fi[iV];
        Fi[iB] = Fi[iV]-Rhs*l0*(l1-l0);
        //Fi[iB] = 0.99999*Fi[iV];
        //Fi[iB] = ((Fi[iV2]-Fi[iV])*l0+Fi[iV]*l2-Fi[iV2]*l1)/dl2;//linear_extrapolation

    if(v==2)///dielectric
        Fi[iB] = (Fi[iV]+Vext*Eps*dl/Dep+4.0*pi*Sig*dl)/(1.0+Eps*dl/Dep);
}

