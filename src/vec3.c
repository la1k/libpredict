#include "vec3.h"
#include <math.h>

void vec3_set(double v[3], double x, double y, double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

double vec3_length(const double v[3])
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double vec3_dot(const double v[3], const double u[3])
{
	return (v[0]*u[0] + v[1]*u[1] + v[2]*u[2]);
}

void vec3_mul_scalar(double v[3], double a)
{
	v[0] *= a;
	v[1] *= a;
	v[2] *= a;
}

void vec3_sub(const double v1[3], const double v2[3], double *r)
{
	r[0] = v1[0] - v2[0];
	r[1] = v1[1] - v2[1];
	r[2] = v1[2] - v2[2];
}
