#include "basic_function.h"

int K_com1 = 0.0;
double r_com1 = 2.0 / (K_com1 + 2) + 1.0;
int K_com2 = 0.0;
double r_com2 = 2.0 / (K_com2 + 2) + 1.0;
double c1_euler = 0.0;
double c2_euler = 0.0;
TAU_TYPE tau_type=Euler; // viscous problem or non-viscous problem

// to prepare the basic element for moment calculation
MMDF1d::MMDF1d() { u = 1.0; lambda = 1.0; flag = 1;}
MMDF1d::MMDF1d(double u_in, double lambda_in, int flag_in)
{
	u = u_in;
	lambda = lambda_in;
	flag = flag_in;
	calcualte_MMDF1d();
}
void MMDF1d::calcualte_MMDF1d()
{
	uwhole[0] = 1;
	uwhole[1] = u;
	uplus[0] = 0.5 * Alpha(lambda, -u);
	uminus[0] = 0.5 * Alpha(lambda, u);
	uplus[1] = u * uplus[0] + 0.5 * Beta(lambda, u);
	uminus[1] = u * uminus[0] - 0.5 * Beta(lambda, u);
	for (int i = 2; i <= 9; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) / lambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 9; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) / lambda * uplus[i - 2];
		uminus[i] = u * uminus[i - 1] + 0.5 * (i - 1) / lambda * uminus[i - 2];
	}

    int K;
    if (flag == 1) {K = K_com1;}
    else {K = K_com2;}

	xi2 = 0.5 * K / lambda;
	xi4 = 0.25 * (K * K + 2 * K) / (lambda * lambda);
	//xi6?? how to calculate
	xi6 = 0.5 * (K + 4) / lambda * xi4;

	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					upxi[i][k] = uplus[i];
					unxi[i][k] = uminus[i];
				}
				if (k == 1)
				{
					upxi[i][k] = uplus[i] * xi2;
					unxi[i][k] = uminus[i] * xi2;
				}
				if (k == 2)
				{
					upxi[i][k] = uplus[i] * xi4;
					unxi[i][k] = uminus[i] * xi4;
				}
				if (k == 3)
				{
					upxi[i][k] = uplus[i] * xi6;
					unxi[i][k] = uminus[i] * xi6;
				}
			}
		}
	}
	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					uxi[i][k] = uwhole[i];
				}
				if (k == 1)
				{
					uxi[i][k] = uwhole[i] * xi2;
				}
				if (k == 2)
				{
					uxi[i][k] = uwhole[i] * xi4;
				}
				if (k == 3)
				{
					uxi[i][k] = uwhole[i] * xi6;
				}
			}
		}
	}
}

double Alpha(double lambda, double u)
{
    return erfc(sqrt(lambda) * u);
}

double Beta(double lambda, double u)
{
	double pi = 3.14159265358979323846;
    return exp(-lambda * u * u) / sqrt(pi * lambda);
}

double Lambda(double density, double u, double densityE, int flag)
{
    int K;
    if (flag == 1){K = K_com1;}
    else{K = K_com2;}

    if (density == 0)
    {
        return 0;
    }
    else
    {
        return (K + 1.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u)));
    }
}

double U(double density, double q_densityu)
{
    if (density == 0) {return 0;}
    else {return q_densityu / density;}
}

void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.uxi[no_u][no_xi] + a[1] * m.uxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.uxi[no_u + 1][no_xi] + a[1] * m.uxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]) +
		a[1] * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.uxi[no_u + 4][no_xi] + m.uxi[no_u][no_xi + 2] + 2 * m.uxi[no_u + 2][no_xi + 1]));

}

void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.upxi[no_u][no_xi] + a[1] * m.upxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.upxi[no_u + 1][no_xi] + a[1] * m.upxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]) +
		a[1] * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.upxi[no_u + 4][no_xi] + m.upxi[no_u][no_xi + 2] + 2 * m.upxi[no_u + 2][no_xi + 1]));
}

void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.unxi[no_u][no_xi] + a[1] * m.unxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.unxi[no_u + 1][no_xi] + a[1] * m.unxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]) +
		a[1] * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.unxi[no_u + 4][no_xi] + m.unxi[no_u][no_xi + 2] + 2 * m.unxi[no_u + 2][no_xi + 1]));
}

//solution of matrix equation b=Ma 1D
void Microslope(Point1d& com1, Point1d& com2)
{
    // prim[] for rho, U, and lambda
    double rho1 = com1.Ulambda[0], rho2 = com2.Ulambda[0];
    double U = com1.Ulambda[1], lambda = com1.Ulambda[2];
    double ome1, ome2, ome3, ome4;
    ome1 = com1.der1[0];
    ome2 = com2.der1[0];
    ome3 = com1.der1[1] + com2.der1[1];
    ome4 = com1.der1[2] + com2.der1[2];

    double psi1, psi2;
    psi1 = ome3 - U * (ome1 + ome2);
    psi2 = 2 * ome4 - (U * U + (K_com1 + 1.0) / (2 * lambda)) * ome1
         - (U * U + (K_com2 + 1.0) / (2 * lambda)) * ome2;

    double p, n;
    p = (2 * lambda * lambda * (psi2 - 2 * U * psi1))
        / ((K_com1 + 1.0) * rho1 + (K_com2 + 1.0) * rho2);
    com1.ax[2] = 2.0 * p;
    com2.ax[2] = 2.0 * p;

    n = (2.0 * lambda / (rho1 + rho2)) * (psi1 - (rho1 + rho2) * U * p / lambda);
    com1.ax[1] = n;
    com2.ax[1] = n;

    com1.ax[0] = (ome1 - rho1 * U * n - rho1 * (U * U + (K_com1 + 1.0) / (2.0 * lambda)) * p) / rho1;
    if (rho1 == 0) { com1.ax[0] = 0; }
    com2.ax[0] = (ome2 - rho2 * U * n - rho2 * (U * U + (K_com2 + 1.0) / (2.0 * lambda)) * p) / rho2;
    if (rho2 == 0) { com2.ax[0] = 0; }
}

Flux1d** Setflux_array(Block1d block)
{
	Flux1d** fluxes = new Flux1d * [block.nx + 1];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i <= block.nx; i++)
	{
		// for m th step time marching schemes, m subflux needed
		fluxes[i] = new Flux1d[block.stages];
		for (int j = 0; j < block.stages; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				fluxes[i][j].f[k] = 0.0;
				fluxes[i][j].derf[k] = 0.0;
			}
		}
	}
	if (fluxes == 0)
	{
		cout << "memory allocation failed for muli-stage flux";
		return NULL;
	}
	cout << "the memory for muli-stage flux has been allocated..." << endl;
	return fluxes;
}

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes1, Flux1d** fluxes2)
{
	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		fluids[i].dx = block.dx; //cell size
		fluids[i].cx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
	}
	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		interfaces[i].x = block.left + (i - block.ghost) * block.dx;
		interfaces[i].comp1.left.x = interfaces[i].x;
		interfaces[i].comp1.right.x = interfaces[i].x;
		interfaces[i].comp1.center.x = interfaces[i].x;
		interfaces[i].flux_comp1 = fluxes1[i];

        interfaces[i].comp2.left.x = interfaces[i].x;
		interfaces[i].comp2.right.x = interfaces[i].x;
		interfaces[i].comp2.center.x = interfaces[i].x;
		interfaces[i].flux_comp2 = fluxes2[i];
	}
}

double U(double rhoU1, double rhoU2, double rho1, double rho2)
{
    return (rhoU1 + rhoU2) / (rho1 + rho2);
}

double Lambda(double rho1, double rho2, double E1, double E2, double U)
{
    return (1.0 / 4.0) * ((K_com1 + 1) * rho1 + (K_com2 + 1) * rho2)
               / (E1 + E2- 0.5 * (rho1 + rho2) * U * U);
}

double Pressure(double density, double densityu, double densityE, double r)
{
    if (density == 0)
    {
        return 0;
    }
    else
    {
        return (r - 1.0) * (densityE - 0.5 * densityu * densityu / density);
    }
}

void Share_speed_and_temperature(Fluid1d* fluids, Block1d block)
{
    for (int i = 0; i < block.nx; i++)
	{
		double velocity, lambda, rho_com1, rho_com2;
        rho_com1 = fluids[i].comp1.convar[0];
        rho_com2 = fluids[i].comp2.convar[0];
        velocity = U(fluids[i].comp1.convar[1], fluids[i].comp2.convar[1], rho_com1, rho_com2);
        lambda = Lambda(rho_com1, rho_com2, fluids[i].comp1.convar[2], fluids[i].comp2.convar[2], velocity);
        fluids[i].comp1.convar[1] = rho_com1 * velocity;
        fluids[i].comp2.convar[1] = rho_com2 * velocity;
        fluids[i].comp1.convar[2] = 0.5 * rho_com1 * (velocity * velocity + (K_com1 + 1.0) / (2.0 * lambda));
        fluids[i].comp2.convar[2] = 0.5 * rho_com2 * (velocity * velocity + (K_com2 + 1.0) / (2.0 * lambda));
	}
}

void Ulambda1d(Interface1d* interfaces, Block1d block)
{
    for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		double U_left, U_right, lambda_left, lambda_right;
        double r_left1, r_left2, r_right1, r_right2;
        r_left1 = interfaces[i].comp1.left.convar[0];
        r_left2 = interfaces[i].comp2.left.convar[0];
        r_right1 = interfaces[i].comp1.right.convar[0];
        r_right2 = interfaces[i].comp2.right.convar[0];

        interfaces[i].comp1.left.Ulambda[0]  = interfaces[i].comp1.left.convar[0];
        interfaces[i].comp2.left.Ulambda[0]  = interfaces[i].comp2.left.convar[0];
        interfaces[i].comp1.right.Ulambda[0] = interfaces[i].comp1.right.convar[0];
        interfaces[i].comp2.right.Ulambda[0] = interfaces[i].comp2.right.convar[0];

        U_left  = (interfaces[i].comp1.left.convar[1] + interfaces[i].comp2.left.convar[1])
                / (r_left1 + r_left2);
        U_right = (interfaces[i].comp1.right.convar[1] + interfaces[i].comp2.right.convar[1])
                / (r_right1 + r_right2);
        interfaces[i].comp1.left.Ulambda[1]  = U_left;
        interfaces[i].comp2.left.Ulambda[1]  = U_left;
        interfaces[i].comp1.right.Ulambda[1] = U_right;
        interfaces[i].comp2.right.Ulambda[1] = U_right;
        // same velocity 
        interfaces[i].comp1.left.convar[1]  = U_left * r_left1;
        interfaces[i].comp2.left.convar[1]  = U_left * r_left2;
        interfaces[i].comp1.right.convar[1] = U_right * r_right1;
        interfaces[i].comp2.right.convar[1] = U_right * r_right2;

        lambda_left  = Lambda(r_left1, r_left2,
                    interfaces[i].comp1.left.convar[2], interfaces[i].comp2.left.convar[2], U_left);
        lambda_right = Lambda(r_right1, r_right2,
                    interfaces[i].comp1.right.convar[2], interfaces[i].comp2.right.convar[2], U_right);
        interfaces[i].comp1.left.Ulambda[2]  = lambda_left;
        interfaces[i].comp2.left.Ulambda[2]  = lambda_left;
        interfaces[i].comp1.right.Ulambda[2] = lambda_right;
        interfaces[i].comp2.right.Ulambda[2] = lambda_right;
        // same temperature
        interfaces[i].comp1.left.convar[2]  = 0.5 * r_left1 * (U_left * U_left + (K_com1 + 1.0) / (2.0 * lambda_left));
        interfaces[i].comp2.left.convar[2]  = 0.5 * r_left2 * (U_left * U_left + (K_com2 + 1.0) / (2.0 * lambda_left));
        interfaces[i].comp1.right.convar[2] = 0.5 * r_right1 * (U_right * U_right + (K_com1 + 1.0) / (2.0 * lambda_right));
        interfaces[i].comp2.right.convar[2] = 0.5 * r_right2 * (U_right * U_right + (K_com2 + 1.0) / (2.0 * lambda_right));
    }
}

void Ulambda_center(Point1d& center1, Point1d& center2)
{
    double U, lambda;
    double r_center1, r_center2;
    r_center1 = center1.convar[0];
    r_center2 = center2.convar[0];

    center1.Ulambda[0]  = r_center1;
    center2.Ulambda[0]  = r_center2;

    U  = (center1.convar[1] + center2.convar[1]) / (r_center1 + r_center2);
    center1.Ulambda[1]  = U;
    center2.Ulambda[1]  = U;
    // same velocity 
    center1.convar[1]  = U * r_center1;
    center2.convar[1]  = U * r_center2;

    lambda  = Lambda(r_center1, r_center2,
                center1.convar[2], center2.convar[2], U);
    center1.Ulambda[2]  = lambda;
    center2.Ulambda[2]  = lambda;
    // same temperature
    center1.convar[2]  = 0.5 * r_center1 * (U * U + (K_com1 + 1.0) / (2.0 * lambda));
    center2.convar[2]  = 0.5 * r_center2 * (U * U + (K_com2 + 1.0) / (2.0 * lambda));
}

double Get_Tau_NS(double density0, double lambda0)
{
	if (tau_type == Euler)
	{
		return 0.0;
	}
	// else
	// {
	// 	// NS
	// 	if (Mu > 0.0)
	// 	{
	// 		//cout << "here" << endl;
	// 		return 2.0 * Mu * lambda0 / density0;
	// 	}
	// 	else if (Nu > 0.0)
	// 	{
	// 		return 2 * Nu * lambda0;
	// 	}
	// 	else
	// 	{
	// 		return 0.0;
	// 	}
	// }
}

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt)
{
	if (tau_type == Euler)
	{
		if (c1_euler <= 0 && c2_euler <= 0)
		{
			return 0.0;
		}
		else
		{
			double C = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
			return c1_euler * dt + dt * C;
		}
	}
	else if (tau_type == NS)
	{
		double tau_n = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right) * dt;
		if (tau_n != tau_n)
		{
			tau_n = 0.0;
		}
		// if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		// {
		// 	return 0.0;
		// }
		// else
		// {
		// 	if (Mu > 0.0)
		// 	{
		// 		return tau_n + 2.0 * Mu * lambda0 / density0;
		// 	}
		// 	else if (Nu > 0.0)
		// 	{
		// 		return tau_n + 2 * Nu * lambda0;
		// 	}
		// 	else
		// 	{
		// 		return 0.0;
		// 	}
		// }
	}
	else
	{
		return 0.0;
	}
}

void Viscous_tau(Point1d& left1, Point1d& right1, Point1d& center1, Point1d& left2, Point1d& right2, Point1d& center2, double dt)
{
    double rho0 = center1.Ulambda[0] + center2.Ulambda[0];
    double lambda0 = center1.Ulambda[2];
    center1.tau = Get_Tau_NS(rho0, lambda0);
    center2.tau = center1.tau;

    double rhol = left1.Ulambda[0] + left2.Ulambda[0];
    double rhor = right1.Ulambda[0] + right2.Ulambda[0];
    double lambdal = left1.Ulambda[2];
    double lambdar = right1.Ulambda[2];
    center1.tau_num = Get_Tau(rhol, rhor, rho0, lambdal, lambdar, lambda0, dt);
    center2.tau_num = center1.tau_num;
}

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block)
{
    Share_speed_and_temperature(fluids, block);
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fluids[i].comp1.convar_old[j] = fluids[i].comp1.convar[j];
            fluids[i].comp2.convar_old[j] = fluids[i].comp2.convar[j];
		}
	}
}

void Convar_to_ULambda_1d(double* primvar, double convar[3], int flag)
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Lambda(convar[0], primvar[1], convar[2], flag);
}