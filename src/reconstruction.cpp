#include "reconstruction.h"

// one-dimensional problem
Reconstruction_within_Cell cellreconstruction = Vanleer; //initialization
Reconstruction_forG0 g0reconstruction = Center_collision; //initialization
Reconstruction_variable reconstruction_variable = conservative; //initialization
WENOtype wenotype = wenojs; //initialization
double df_thres = 3.0;
bool is_reduce_order_warning = false; //initialization
bool is_using_df_factor = false;
bool smooth = false;

void Reconstruction_within_cell(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
        // component 1
		cellreconstruction(interfaces[i].comp1.right, interfaces[i + 1].comp1.left, &fluids[i], block, 1);
		// component 2
		cellreconstruction(interfaces[i].comp2.right, interfaces[i + 1].comp2.left, &fluids[i], block, 2);
		// check if the lambda <= 0 or not
        Check_Order_Reduce(interfaces[i].comp1.right, interfaces[i + 1].comp1.left, interfaces[i].comp2.right, interfaces[i + 1].comp2.left, &fluids[i]);
    }
}

void Check_Order_Reduce(Point1d& left_com1, Point1d& right_com1, Point1d& left_com2, Point1d& right_com2,
                        Fluid1d* fluid)
{
    double lambda_left, lambda_right;
    double U_left = (left_com1.convar[1] + left_com2.convar[1]) / (left_com1.convar[0] + left_com2.convar[0]);
    double U_right = (right_com1.convar[1] + right_com2.convar[1]) / (right_com1.convar[0] + right_com2.convar[0]);

    lambda_left = Lambda(left_com1.convar[0], left_com2.convar[0], left_com1.convar[2], left_com2.convar[2], U_left);
	lambda_right = Lambda(right_com1.convar[0], right_com2.convar[0], right_com1.convar[2], right_com2.convar[2], U_right);

	//if lambda <0, then reduce to the first order
	if (lambda_left <= 0 || lambda_right <= 0)
	{
		for (int m = 0; m < 3; m++)
		{
            // component 1
            left_com1.convar[m] = fluid[0].comp1.convar[m];
            right_com1.convar[m] = fluid[0].comp1.convar[m];
            left_com1.der1[m] = 0.0;
			right_com1.der1[m] = 0.0;
            // component 2
            left_com2.convar[m] = fluid[0].comp2.convar[m];
            right_com2.convar[m] = fluid[0].comp2.convar[m];
            left_com2.der1[m] = 0.0;
			right_com2.der1[m] = 0.0;
		}
	}
}

void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block, int flag)
{
	//Note: function by vanleer reconstruction
	//the final updated variables is conservative variables
	Fluid1d wn1 = fluids[-1];
	Fluid1d w = fluids[0];
	Fluid1d wp1 = fluids[1];
	double splus[3], sminus[3];

	for (int i = 0; i < 3; i++)
	{
        if (flag == 1)
        {
            // component 1
            splus[i] = (wp1.comp1.convar[i] - w.comp1.convar[i]) / ((wp1.dx + w.dx) / 2.0);
		    sminus[i] = (w.comp1.convar[i] - wn1.comp1.convar[i]) / ((wn1.dx + w.dx) / 2.0);
        }
        else
        {
            // component 2
            splus[i] = (wp1.comp2.convar[i] - w.comp2.convar[i]) / ((wp1.dx + w.dx) / 2.0);
		    sminus[i] = (w.comp2.convar[i] - wn1.comp2.convar[i]) / ((wn1.dx + w.dx) / 2.0);
        }

		if ((splus[i] * sminus[i]) > 0)
		{
			left.der1[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
			right.der1[i] = left.der1[i];
		}
		else
		{
			left.der1[i] = 0.0;
			right.der1[i] = 0.0;
		}
		if (flag == 1)
        {
            // component 1
            left.convar[i] = w.comp1.convar[i] - 0.5 * w.dx * left.der1[i];
		    right.convar[i] = w.comp1.convar[i] + 0.5 * w.dx * right.der1[i];
        }
        else
        {
            // component 2
            left.convar[i] = w.comp2.convar[i] - 0.5 * w.dx * left.der1[i];
		    right.convar[i] = w.comp2.convar[i] + 0.5 * w.dx * right.der1[i];
        }
	}
	//Check_Order_Reduce(left, right, fluids[0], flag);
}

// interface center reconstruction (only 2nd-order GKS is needed)
void Reconstruction_forg0(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; ++i)
	{
		// left coefficient Ma = b
        Microslope(interfaces[i].comp1.left, interfaces[i].comp2.left);
		// right coefficient Ma = b
        Microslope(interfaces[i].comp1.right, interfaces[i].comp2.right);
		// component 1
		g0reconstruction(interfaces[i].comp1.left, interfaces[i].comp1.right, interfaces[i].comp1.center, &fluids[i - 1], 1);
		// component 2
		g0reconstruction(interfaces[i].comp2.left, interfaces[i].comp2.right, interfaces[i].comp2.center, &fluids[i - 1], 2);
		Ulambda_center(interfaces[i].comp1.center, interfaces[i].comp2.center);
		// center coefficient Ma = b
        Microslope(interfaces[i].comp1.center, interfaces[i].comp2.center);
    }
}

void Center_collision(Point1d& left, Point1d& right, Point1d& center, Fluid1d* fluids, int flag)
{
	double convar_left[3], convar_right[3];
    double prim_left[3], prim_right[3]; //rho, U, lambda
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = left.convar[i];
		convar_right[i] = right.convar[i];
        prim_left[i] = left.Ulambda[i];
        prim_right[i] = right.Ulambda[i];
	}

	MMDF1d ml(prim_left[1], prim_left[2], flag);
	MMDF1d mr(prim_right[1], prim_right[2], flag);
	
	double unit[3]{ 1.0, 0.0, 0.0 };

	double gl[3], gr[3];
	GL(0, 0, gl, unit, ml); // gl, means Wl(u>0), by input uint
	GR(0, 0, gr, unit, mr); // gr, means Wr(u<0), by input uint
	double axl[3], axr[3];
    for (int i = 0; i < 3; i++)
    {
        axl[i] = left.ax[i];
        axr[i] = right.ax[i];
    }
	// Microslope(axl, left.der1, prim_left, flag); // axl, means a coefficient indicating slope
	// Microslope(axr, right.der1, prim_right, flag); // axr, means a coefficient indicating slope
	double ax0l[3], ax0r[3];
	GL(0, 0, ax0l, axl, ml);  // ax0l, means Wlx(u>0), by input axl
	GR(0, 0, ax0r, axr, mr); // ax0r, means Wrx(u<0), by input axr
	for (int i = 0; i < 3; i++)
	{
		center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
		center.der1[i] = convar_left[0] * ax0l[i] + convar_right[0] * ax0r[i];
	}
}
