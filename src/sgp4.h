#ifndef SGP4_H_
#define SGP4_H_

#include "unsorted.h"
#include "vec3.h"

struct _sgp4 {
	///Phase?
	double phase;
	///Initialized flag:
	int initialized;
	///Simple flag
	int simpleFlag;
	///Static variables from original SGP4()
	double aodp, aycof, c1, c4, c5, cosio, d2, d3, d4, delmo,
	omgcof, eta, omgdot, sinio, xnodp, sinmo, t2cof, t3cof, t4cof,
	t5cof, x1mth2, x3thm1, x7thm1, xmcof, xmdot, xnodcf, xnodot, xlcof;
	
};
	
void sgp4_init(struct _sgp4 *m);
void sgp4_predict(struct _sgp4 *m, double tsince, tle_t *tle, double pos[3], double vel[3]);

#endif
