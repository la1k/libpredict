#ifndef SGP4_H_
#define SGP4_H_

#include <predict/predict.h>

/**
 * Parameters relevant for SGP4 (simplified general perturbations) orbital model.
 **/
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

/**
 * Initialize SGP4 model parameters.  
 *
 * \param m Struct to initialize
 **/
void sgp4_init(struct _sgp4 *m);

/**
 * Predict ECI position and velocity of near-earth orbit (period < 225 minutes) according to SGP4 model and the given orbital parameters. 
 *
 * \param m SGP4 model parameters
 * \param tsince Time since epoch of TLE in minutes
 * \param tle TLE parameters
 * \param pos Output ECI position in meters
 * \param vel Output velocity in m/s
 * \copyright GPLv2+
 **/
void sgp4_predict(struct _sgp4 *m, double tsince, predict_tle_t *tle, double pos[3], double vel[3]);

#endif
