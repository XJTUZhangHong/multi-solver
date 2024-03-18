#pragma once
#include "fluid_mesh.h"
#include <assert.h>

extern int K_com1;
extern double r_com1;
extern int K_com2;
extern double r_com2;
extern double c1_euler;
extern double c2_euler;
enum TAU_TYPE { Euler, NS};
extern TAU_TYPE tau_type;

//basic gks function
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;
    int flag;
public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in, int flag);
	void calcualte_MMDF1d();
};

double Alpha(double lambda, double u);

double Beta(double lambda, double u);

double U(double density, double q_densityu);

double Lambda(double density, double u, double densityE, int flag);

void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void Microslope(Point1d& com1, Point1d& com2);

Flux1d** Setflux_array(Block1d block);

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes1, Flux1d** fluxes2);

double U(double rhoU1, double rhoU2, double rho1, double rho2);

double Lambda(double rho1, double rho2, double E1, double E2, double U);

double Pressure(double density, double densityu, double densityE, double r);

void Share_speed_and_temperature(Fluid1d* fluids, Block1d block);

void Ulambda1d(Interface1d* interfaces, Block1d block);

void Ulambda_center(Point1d& center1, Point1d& center2);

double Get_Tau_NS(double density0, double lambda0);

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt);

void Viscous_tau(Point1d& left1, Point1d& right1, Point1d& center1, Point1d& left2, Point1d& right2, Point1d& center2, double dt);

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block);

void Convar_to_ULambda_1d(double* primvar, double convar[3], int flag);
