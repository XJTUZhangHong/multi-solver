#include "flux_function.h"

// one-dimensional problem
GKS1d_type gks1dsolver = nothing; //initialization
Flux_function flux_function = GKS; //initialization

void Calculate_flux(Flux1d** fluxes1, Flux1d** fluxes2, Interface1d* interfaces, Block1d& block, int stage)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
	{
        // component 1
		flux_function(fluxes1[i][stage], interfaces[i].comp1.left, interfaces[i].comp1.right, interfaces[i].comp1.center, block.dt, 1);
		// component 2
        flux_function(fluxes2[i][stage], interfaces[i].comp2.left, interfaces[i].comp2.right, interfaces[i].comp2.center, block.dt, 2);
	}
}

void GKS(Flux1d& flux, Point1d& left, Point1d& right, Point1d& center, double dt, int flag)
{
	if (gks1dsolver == nothing)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	double Flux[2][3];
	//change conservative variables to rho u lambda
	double convar_left[3], convar_right[3], convar0[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = left.convar[i];
		convar_right[i] = right.convar[i];
		convar0[i] = center.convar[i];
	}

	double prim_left[3], prim_right[3], prim0[3];
	for (int i = 0; i < 3; i++)
    {
        prim_left[i] = left.Ulambda[i];
        prim_right[i] = right.Ulambda[i];
        prim0[i] = center.Ulambda[i];
    }
	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[2]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm (f0)
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd (g0)
	t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part

	if (gks1dsolver == kfvs1st)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks1dsolver == kfvs2nd)
	{
		t[0] = dt;
		t[1] = -dt * dt / 2.0;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}
	MMDF1d ml(prim_left[1], prim_left[2], flag);
	MMDF1d mr(prim_right[1], prim_right[2], flag);

	double unit[3] = { 1, 0.0, 0.0 };

	double glu[3], gru[3];
	GL(1, 0, glu, unit, ml);
	GR(1, 0, gru, unit, mr);

	//only one part, the kfvs1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}

	if (gks1dsolver == kfvs1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF1d m0(prim0[1], prim0[2], flag);

	double g0u[3];
	G(1, 0, g0u, unit, m0);

	//the equ g0u part, the gks1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
	}

	if (gks1dsolver == gks1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1left[3], der1right[3];
	for (int i = 0; i < 3; i++)
	{
		der1left[i] = left.der1[i];
		der1right[i] = right.der1[i];
	}

	double alx[3];
    for (int i = 0; i < 3; i++)
    {
        alx[i] = left.ax[i];
    }
	// Microslope(alx, der1left, prim_left);

	double alxuul[3];
	GL(2, 0, alxuul, alx, ml);

	double arx[3];
    for (int i = 0; i < 3; i++)
    {
        arx[i] = right.ax[i];
    }
	// Microslope(arx, der1right, prim_right);
	double arxuur[3];
	GR(2, 0, arxuur, arx, mr);

	for (int i = 0; i < 3; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
	}
	if (gks1dsolver == kfvs2nd)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd
	//for t[2] Aru,Alu part
	double alxu[3];
	double arxu[3];

	//take <u> moment for al, ar
	G(1, 0, alxu, alx, ml);
	G(1, 0, arxu, arx, mr);

	double Al[3], Ar[3];
	double der_AL[3], der_AR[3];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 3; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i]);
		der_AR[i] = -prim_right[0] * (arxu[i]);
	}
	// solve the coefficient martix b=ma
    //-----------------------------------//
	// Microslope(Al, der_AL, prim_left);
	// Microslope(Ar, der_AR, prim_right);
    //-----------------------------------//

	//to obtain the Alu and Aru
	double Alul[3];
	double Arur[3];
	GL(1, 0, Alul, Al, ml);
	GR(1, 0, Arur, Ar, mr);

	for (int i = 0; i < 3; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
	}

	// for t[4] a0xuu part

	double a0x[3];
	double der1[3];

	for (int i = 0; i < 3; i++)
	{
		der1[i] = center.der1[i];
        a0x[i] = center.ax[i];
	}

	//solve the microslope
	//Microslope(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[3];
	G(1, 0, a0xu, a0x, m0);
	//a0x <u^2> moment
	double a0xuu[3];
	G(2, 0, a0xuu, a0x, m0);

	for (int i = 0; i < 3; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i]);
	}


	// for t[5] A0u part
	double derA0[3];

	for (int i = 0; i < 3; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i]);
	}
	double A0[3];
    //-----------------------------------//
	// Microslope(A0, derA0, prim0);
    //-----------------------------------//
	double A0u[3];
	G(1, 0, A0u, A0, m0);
	for (int i = 0; i < 3; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (gks1dsolver == gks2nd && timecoe_list == S1O1)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	if (gks1dsolver == gks2nd)
	{
		double dt2 = 0.5 * dt;
		tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num * (1 - eta); // this refers glu, gru part
		t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
		t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

		for (int i = 0; i < 3; i++)
		{
			// t0 part
			Flux[1][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
			// t1 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
			// t2 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}

		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = (4.0 * Flux[1][i] - Flux[0][i]);
			flux.derf[i] = 4.0 * (Flux[0][i] - 2.0 * Flux[1][i]);
		}

		return;
	}
	else
	{
		cout << "no valid solver specify" << endl;
		exit(0);
	}
}
