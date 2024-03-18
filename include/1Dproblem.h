#pragma once
#include "flux_function.h"
#include "output.h"
#include "reconstruction.h"
#include "time_advance.h"
#include "boundary_condition.h"

void SodTubeProblem();

void ICforSodTube(Fluid1d* fluids, Block1d block);
