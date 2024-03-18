#pragma once
#include "flux_function.h"

enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
enum WENOtype { linear, wenojs, wenoz };
extern WENOtype wenotype;
extern double df_thres;
extern bool is_reduce_order_warning;
extern bool is_using_df_factor;
extern bool smooth;

// one-dimensional problem
// interface left & right reconstruction
void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void(*Reconstruction_within_Cell)(Point1d &left, Point1d &right, Fluid1d *fluids, Block1d block, int flag);
extern Reconstruction_within_Cell cellreconstruction;
void Check_Order_Reduce(Point1d& left_com1, Point1d& right_com1, Point1d& left_com2, Point1d& right_com2,
                        Fluid1d* fluid);
void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block, int flag);

// interface center reconstruction
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Point1d& left, Point1d& right, Point1d& center, Fluid1d* fluids, int flag);
extern Reconstruction_forG0 g0reconstruction;
void Center_collision(Point1d& left, Point1d& right, Point1d& center, Fluid1d* fluids, int flag);