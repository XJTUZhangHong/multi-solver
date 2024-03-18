#pragma once
#include "time_advance.h"
#include "reconstruction.h"

// one-dimensional problem
enum GKS1d_type{nothing, kfvs1st, kfvs2nd, gks1st, gks2nd, pp_gks}; // positive-preserving gks
extern GKS1d_type gks1dsolver;

void Calculate_flux(Flux1d** fluxes1, Flux1d** fluxes2, Interface1d* interfaces, Block1d& block, int stage);
typedef void(*Flux_function)(Flux1d& flux, Point1d& left, Point1d& right, Point1d& center, double dt, int flag);
extern Flux_function flux_function;
void GKS(Flux1d& flux, Point1d& left, Point1d& right, Point1d& center, double dt, int flag);
