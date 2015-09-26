#ifndef _SDP4_H_
#define _SDP4_H_

#include <predict/predict.h>

/**
 * Parameters for deep space perturbations?
 **/
typedef struct	{
		   	   /* Used by dpinit part of Deep() */
		   double  eosq, sinio, cosio, betao, aodp, theta2,
			   sing, cosg, betao2, xmdot, omgdot, xnodot, xnodp;
	   
			   /* Used by dpsec and dpper parts of Deep() */
		   double  xll, omgadf, xnode, em, xinc, xn, t;
    
		 	   /* Used by thetg and Deep() */
		   double  ds50;
		}  deep_arg_t;

/**
 * Parameters relevant for SDP4 (simplified deep space perturbations) orbital model.
 **/
struct _sdp4 {
	
	//Phase?
	double phase;

	///Have the SDP4 function been initialized?
	int initialized;
	///Lunar terms done?
	int lunarTermsDone;
	///Resonance flag:
	int resonanceFlag;
	///Synchronous flag:
	int synchronousFlag;
	///Do loop flag:
	int loopFlag;
	///Epoch restart flag:
	int epochRestartFlag;

	
	///Static variables from SDP4():
	double x3thm1, c1, x1mth2, c4, xnodcf, t2cof, xlcof,
	aycof, x7thm1;
	deep_arg_t deep_arg;
	
	///Static variables from Deep():
	double thgr, xnq, xqncl, omegaq, zmol, zmos, savtsn, ee2, e3,
	xi2, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2, xh3, sse, ssi, ssg, xi3,
	se2, si2, sl2, sgh2, sh2, se3, si3, sl3, sgh3, sh3, sl4, sgh4, ssl,
	ssh, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, del1,
	del2, del3, fasx2, fasx4, fasx6, xlamo, xfact, xni, atime, stepp,
	stepn, step2, preep, pl, sghs, xli, d2201, d2211, sghl, sh1, pinc,
	pe, shs, zsingl, zcosgl, zsinhl, zcoshl, zsinil, zcosil;

	//Variables that are used locally in SDP4(), but also are used to calculate squint angle elsewhere
	double xnodek, xinck;
};

/**
 * Initialize SDP4 model parameters. 
 *
 * \param m Struct to initialize
 **/
void sdp4_init(struct _sdp4 *m);

/** 
 * Predict ECI position and velocity of deep-space orbit (period > 225 minutes) according to SDP4 model and the given orbital parameters. 
 *
 * \param m SDP4 model parameters
 * \param tsince Time since epoch of TLE in minutes
 * \param tle TLE parameters
 * \param pos Output ECI position in meters
 * \param vel Output velocity in m/s
 * \copyright GPLv2+
 **/
void sdp4_predict(struct _sdp4 *m, double tsince, predict_tle_t * tle, double pos[3], double vel[3]);

/**
 * Deep space perturbations. Original Deep() function.
 *
 * \param m SDP4 model parameters
 * \param ientry Behavior flag. 0: Deep space initialization. 1: Deep space secular effects. 2: lunar-solar periodics
 * \param tle TLE parameters
 * \param deep_arg Deep perturbation parameters
 * \copyright GPLv2+
 **/
void sdp4_deep(struct _sdp4 *m, int ientry, predict_tle_t * tle, deep_arg_t * deep_arg);


#endif // ifndef _SDP4_H_
