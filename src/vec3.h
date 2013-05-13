#ifndef VEC3_H_
#define VEC3_H_

void vec3_set(double v[3], double x, double y, double z);
double vec3_length(const double v[3]);
void vec3_mul_scalar(double v[3], double a);
void vec3_sub(const double v1[3], const double v2[3], double *r);
double vec3_dot(const double v[3], const double u[3]);

#endif
