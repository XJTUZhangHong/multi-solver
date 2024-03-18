#include "1Dproblem.h"

void SodTubeProblem()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 100;
	block.ghost = 3;

	double tstop = 0.2;
	block.CFL = 0.5;
	Fluid1d* bcvalue = new Fluid1d[2];
	r_com1 = 1.4;
    r_com2 = 1.2;
    K_com1 = (5.0 - 3.0 * r_com1) / (r_com1 - 1.0) + 2.0;
    K_com2 = (5.0 - 3.0 * r_com2) / (r_com2 - 1.0) + 2.0;

	gks1dsolver = kfvs2nd;
	//g0type = collisionn;
	tau_type = Euler;
	//Smooth = false;
	c1_euler = 0.05;
	c2_euler = 1.0;

	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	//prepare the reconstruction function

	cellreconstruction = Vanleer;
	wenotype = wenoz;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_collision;
	is_reduce_order_warning = true;
	//prepare the flux function
	flux_function = GKS;
	//prepare time marching stratedgy
	timecoe_list = S1O2;
	Initial_stages(block);

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes_com1 = Setflux_array(block);
    Flux1d** fluxes_com2 = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = 0.0;
	block.right = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes_com1, fluxes_com2);
	//ended mesh part

	ICforSodTube(fluids, block);
	//initializing part end

	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	runtime.finish_initial = clock();
	while (block.t < tstop)
	{
		// assume you are using command window,
		// you can specify a running step conveniently
		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output1d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

            // here need to calculate the same velocity and temperature
            Share_speed_and_temperature(fluids, block);
			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);

            // here need to calculate the same velocity and temperature for the reconstructed values
            Ulambda1d(interfaces, block);

			Reconstruction_forg0(interfaces, fluids, block);
			
			//then is solver part
			Calculate_flux(fluxes_com1, fluxes_com2, interfaces, block, i);
			//then is update flux part
			Update(fluids, fluxes_com1, fluxes_com2, block, i);
		}
		// update the compression factor
		//for (int j = 3; j < 402; j++) { cout << fluids[j].convar[0] << endl; }
		block.step++;
		block.t = block.t + block.dt;
		if (block.step % 1000 == 0)
		{
			cout << "step 1000 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) > 0)
		{
			runtime.finish_compute = clock();
			cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
			cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
			output1d(fluids, block);
			//output1d_checking(fluids, interfaces, fluxes, block);
		}
	}
}

void ICforSodTube(Fluid1d* fluids, Block1d block)
{
	for (int i = block.ghost; i <= block.nodex + 2 * block.ghost - 1; i++)
	{
		if ((double)(i - block.ghost + 1) * block.dx <= 0.5)
		{
            // component 1
			fluids[i].comp1.convar[0] = 1.0;
			fluids[i].comp1.convar[1] = 0.0;
			fluids[i].comp1.convar[2] = 2.5;
            // component 2
			fluids[i].comp2.convar[0] = 0.0;
			fluids[i].comp2.convar[1] = 0.0;
			fluids[i].comp2.convar[2] = 0.0;
		}
		else
		{
            // component 2
			fluids[i].comp1.convar[0] = 0.0;
			fluids[i].comp1.convar[1] = 0.0;
			fluids[i].comp1.convar[2] = 0.0;
            // component 2
			fluids[i].comp2.convar[0] = 0.125;
			fluids[i].comp2.convar[1] = 0.0;
			fluids[i].comp2.convar[2] = 0.5;
		}
	}
}
