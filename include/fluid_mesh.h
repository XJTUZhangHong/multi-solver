#pragma once
#include <omp.h>
#include <cmath>
#include <iostream>
using namespace std;

// one-dimensional problem
class Block1d
{
public:
	//bool uniform;
	int ghost;
	int nodex; // mesh number
	int nx;
	int nodex_begin;
	int nodex_end;
	double dx;
	double left;
	double right;
	int stages;
	double timecoefficient[5][5][2];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
};

class Compent_cell_value1d
{
public:
    double primvar[3];
    double convar[3];
    double convar_old[3];
};

class Fluid1d
{
public:
	Compent_cell_value1d comp1;
    Compent_cell_value1d comp2;
	double cx; //center coordinate in x direction
	double dx; //the mesh size dx
	double alpha = 0.0;
};

// remember the flux in a fixed interface,
// for RK method, we only need f
// for 2nd der method, we need derf
class Flux1d
{
public:
	double F[3]; //total flux in dt time
	double f[3]; // the f0 in t=0
	double derf[4]; // the f_t in t=0
};

class Point1d
{
public:
	double convar[3];
	double convar_old[3]; //this one is for compact point-wise reconstruction
	double der1[3];
	double x; // coordinate of a point
	double Ulambda[3];
	double ax[3]; // coefficient Ma=b
	double tau;
	double tau_num;
};

class Component_interface_value1d
{
public:
    Point1d left;
	Point1d center;
	Point1d right;
};

class Interface1d
{
public:
	Component_interface_value1d comp1;
    Component_interface_value1d comp2;
    Flux1d* flux_comp1;
    Flux1d* flux_comp2;
	double x; // coordinate of the interface, equal to point1d.x
};
