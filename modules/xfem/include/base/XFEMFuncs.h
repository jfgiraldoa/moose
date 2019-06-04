//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>
#include "libmesh/point.h"
#include "EFAPoint.h"

using namespace libMesh;

namespace Xfem
{
static const double tol = 1.0e-10;

void dunavant_rule2(const Real * wts,
                    const Real * a,
                    const Real * b,
                    const unsigned int * permutation_ids,
                    unsigned int n_wts,
                    std::vector<Point> & points,
                    std::vector<Real> & weights);

void stdQuadr2D(unsigned int nen, unsigned int iord, std::vector<std::vector<Real>> & sg2);

void wissmannPoints(unsigned int nqp, std::vector<std::vector<Real>> & wss);

void shapeFunc2D(unsigned int nen,
                 std::vector<Real> & ss,
                 std::vector<Point> & xl,
                 std::vector<std::vector<Real>> & shp,
                 Real & xsj,
                 bool natl_flg);

double r8vec_norm(int n, double a[]);

void r8vec_copy(int n, double a1[], double a2[]);

bool r8vec_eq(int n, double a1[], double a2[]);

double r8vec_dot_product(int n, double a1[], double a2[]);

bool line_exp_is_degenerate_nd(int dim_num, double p1[], double p2[]);

int plane_normal_line_exp_int_3d(
    double pp[3], double normal[3], double p1[3], double p2[3], double pint[3]);

double polyhedron_volume_3d(
    double coord[], int order_max, int face_num, int node[], int node_num, int order[]);

void i4vec_zero(int n, int a[]);

void normalizePoint(Point & p);

void normalizePoint(EFAPoint & p);

double r8_acos(double c);

double angle_rad_3d(double p1[3], double p2[3], double p3[3]);
} // namespace Xfem

