#ifndef _SDP4_H_
#define _SDP4_H_

#include "unsorted.h"

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


};


void sdp4_init(struct _sdp4 *m);

/* This function is used to calculate the position and velocity */
/* of deep-space (period > 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s. */
void sdp4_predict(struct _sdp4 *m, double tsince, tle_t * tle, double pos[3], double vel[3]);

///This is the original Deep() function
void sdp4_deep(struct _sdp4 *m, int ientry, tle_t * tle, deep_arg_t * deep_arg);


#endif // ifndef _SDP4_H_
