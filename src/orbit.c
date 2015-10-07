#define _XOPEN_SOURCE 600
#include <math.h>
#include "defs.h"
#include "unsorted.h"
#include "sdp4.h"
#include "sgp4.h"
#include "sun.h"

bool is_eclipsed(const double pos[3], const double sol[3], double *depth);

/**
 * \brief Allocates memory and prepares internal data.
 *
 * This function is a combination of the original PreCalc() and
 * select_ephemeris() functions.
 **/
predict_orbit_t *predict_create_orbit(char *tle[])
{
	// Allocate memory for new orbit struct
	predict_orbit_t *m = (predict_orbit_t*)malloc(sizeof(predict_orbit_t));
	if (m == NULL) return NULL;

	m->time = nan("");
	m->position[0] = nan("");
	m->position[1] = nan("");
	m->position[2] = nan("");
	m->velocity[0] = nan("");
	m->velocity[1] = nan("");
	m->velocity[2] = nan("");
	m->latitude = nan("");
	m->longitude = nan("");
	m->altitude = nan("");
	m->eclipsed = false;
	m->eclipse_depth = nan("");
	m->footprint = nan("");
	m->phase = nan("");
	m->revolutions = 0;

	//Parse TLE
	double tempnum;
	m->tle.catnum = atol(SubString(tle[0],2,6));
	m->tle.setnum = atol(SubString(tle[0],64,67));
	m->tle.year = atoi(SubString(tle[0],18,19));
	m->tle.refepoch = atof(SubString(tle[0],20,31));
	m->tle.incl = atof(SubString(tle[1],8,15));
	m->tle.raan = atof(SubString(tle[1],17,24));
	m->tle.eccn = 1.0e-07*atof(SubString(tle[1],26,32));
	m->tle.argper = atof(SubString(tle[1],34,41));
	m->tle.meanan = atof(SubString(tle[1],43,50));
	m->tle.meanmo = atof(SubString(tle[1],52,62));
	m->tle.drag  = atof(SubString(tle[0],33,42));
	tempnum=1.0e-5*atof(SubString(tle[0],44,49));
	m->tle.nddot6 = tempnum/pow(10.0,(tle[0][51]-'0'));
	tempnum=1.0e-5*atof(SubString(tle[0],53,58));
	m->tle.bstar = tempnum/pow(10.0,(tle[0][60]-'0'));
	m->tle.orbitnum = atof(SubString(tle[1],63,67));
	
	/* Period > 225 minutes is deep space */
	double ao, xnodp, dd1, dd2, delo, a1, del1, r1;
	double temp = twopi/xmnpda/xmnpda;
	double xno = m->tle.meanmo*temp*xmnpda; //from old TLE struct
	dd1=(xke/xno);
	dd2=tothrd;
	a1=pow(dd1,dd2);
	r1=cos(m->tle.incl*M_PI/180.0);
	dd1=(1.0-m->tle.eccn*m->tle.eccn);
	temp=ck2*1.5f*(r1*r1*3.0-1.0)/pow(dd1,1.5);
	del1=temp/(a1*a1);
	ao=a1*(1.0-del1*(tothrd*.5+del1*(del1*1.654320987654321+1.0)));
	delo=temp/(ao*ao);
	xnodp=xno/(delo+1.0);

		
	/* Select a deep-space/near-earth ephemeris */
	if (twopi/xnodp/xmnpda >= 0.15625) {
		
		m->ephemeris = EPHEMERIS_SDP4;

		// Allocate memory for ephemeris data
		m->ephemeris_data = malloc(sizeof(struct _sdp4));

		if (m->ephemeris_data == NULL) {
			predict_destroy_orbit(m);
			return NULL;
		}
		// Initialize ephemeris data structure
		sdp4_init((struct _sdp4*)m->ephemeris_data);

	}else {
		m->ephemeris = EPHEMERIS_SGP4;
		
		// Allocate memory for ephemeris data
		m->ephemeris_data = malloc(sizeof(struct _sgp4));

		if (m->ephemeris_data == NULL) {
			predict_destroy_orbit(m);
			return NULL;
		}
		// Initialize ephemeris data structure
		sgp4_init((struct _sgp4*)m->ephemeris_data);
	}

	return m;
}

void predict_destroy_orbit(predict_orbit_t *orbit)
{
	if (orbit == NULL) return;

	if (orbit->ephemeris_data != NULL) {
		free(orbit->ephemeris_data);
	}

	free(orbit);
}

/**
 * \brief Is this satellite considered geostationary?
 *
 * This function returns true if the satellite is geostationary.
 **/
bool predict_is_geostationary(const predict_orbit_t *m)
{
	if (fabs(m->tle.meanmo-1.0027)<0.0002) {
		return true;
	}else {
		return false;
	}
}

double predict_apogee(const predict_orbit_t *m)
{
	double sma = 331.25*exp(log(1440.0/m->tle.meanmo)*(2.0/3.0));
	return sma*(1.0+m->tle.eccn)-xkmper;
}
		
double predict_perigee(const predict_orbit_t *m)
{
	double xno = m->tle.meanmo*twopi/xmnpda;
	double a1=pow(xke/xno,tothrd);
	double cosio=cos(m->tle.incl*M_PI/180.0);
	double theta2=cosio*cosio;
	double x3thm1=3*theta2-1.0;
	double eosq=m->tle.eccn*m->tle.eccn;
	double betao2=1.0-eosq;
	double betao=sqrt(betao2);
	double del1=1.5*ck2*x3thm1/(a1*a1*betao*betao2);
	double ao=a1*(1.0-del1*(0.5*tothrd+del1*(1.0+134.0/81.0*del1)));
	double delo=1.5*ck2*x3thm1/(ao*ao*betao*betao2);
	double aodp=ao/(1.0-delo);

	return (aodp*(1-m->tle.eccn)-ae)*xkmper;
}

bool predict_aos_happens(const predict_orbit_t *m, double latitude)
{
	/* This function returns true if the satellite pointed to by
	   "x" can ever rise above the horizon of the ground station. */

	double lin, sma, apogee;

	if (m->tle.meanmo==0.0)
		return false;
	else
	{
		lin = m->tle.incl;

		if (lin >= 90.0) lin = 180.0-lin;

		sma = 331.25*exp(log(1440.0/m->tle.meanmo)*(2.0/3.0));
		apogee = sma*(1.0+m->tle.eccn)-xkmper;

		if ((acos(xkmper/(apogee+xkmper))+(lin*M_PI/180.0)) > fabs(latitude*M_PI/180.0))
			return true;
		else
			return false;
	}
}

/* This is the stuff we need to do repetitively while tracking. */
/* This is the old Calc() function. */
int predict_orbit(predict_orbit_t *m, double utc)
{
	/* Set time to now if now time is provided: */
	if (utc == 0) utc = predict_to_julian(time(NULL));
	
	/* Satellite position and velocity vectors */
	vec3_set(m->position, 0, 0, 0);
	vec3_set(m->velocity, 0, 0, 0);

	m->time = utc;
	double julTime = utc + 2444238.5;

	/* Convert satellite's epoch time to Julian  */
	/* and calculate time since epoch in minutes */
	double epoch = 1000.0*m->tle.year + m->tle.refepoch;
	double jul_epoch = Julian_Date_of_Epoch(epoch);
	double tsince = (julTime - jul_epoch)*xmnpda;

	/* Call NORAD routines according to deep-space flag. */
	switch (m->ephemeris) {
		case EPHEMERIS_SDP4:
			sdp4_predict((struct _sdp4*)m->ephemeris_data, tsince, &m->tle, m->position, m->velocity);
			m->phase = ((struct _sdp4*)m->ephemeris_data)->phase;
			break;
		case EPHEMERIS_SGP4:
			sgp4_predict((struct _sgp4*)m->ephemeris_data, tsince, &m->tle, m->position, m->velocity);
			m->phase = ((struct _sgp4*)m->ephemeris_data)->phase;
			break;
		default:
			//Panic!
			return -1;
	}

	/* TODO: Remove? Scale position and velocity vectors to km and km/sec */
	Convert_Sat_State(m->position, m->velocity);

	/* Calculate satellite Lat North, Lon East and Alt. */
	geodetic_t sat_geodetic;
	Calculate_LatLonAlt(utc, m->position, &sat_geodetic);

	m->latitude = sat_geodetic.lat;
	m->longitude = sat_geodetic.lon;
	m->altitude = sat_geodetic.alt;

	// Calculate solar position
	double solar_vector[3];
	sun_predict(m->time, solar_vector);

	// Find eclipse depth and if sat is eclipsed
	m->eclipsed = is_eclipsed(m->position, solar_vector, &m->eclipse_depth);

	// Calculate footprint
	m->footprint = 2.0*xkmper*acos(xkmper/(xkmper + m->altitude));
	
	// Calculate current number of revolutions around Earth
	double temp = twopi/xmnpda/xmnpda;
	double age = julTime - jul_epoch;
	double xno = m->tle.meanmo*temp*xmnpda;
	double xmo = m->tle.meanan * M_PI / 180.0;
	m->revolutions = (long)floor((xno*xmnpda/(M_PI*2.0) + age*m->tle.bstar)*age + xmo/(2.0*M_PI)) + m->tle.orbitnum;

	return 0;
}

bool predict_decayed(const predict_orbit_t *orbit)
{
	double time = orbit->time;
	double satepoch;
	satepoch=DayNum(1,0,orbit->tle.year)+orbit->tle.refepoch;

	bool has_decayed = false;
	if (satepoch + ((16.666666 - orbit->tle.meanmo)/(10.0*fabs(orbit->tle.drag))) < time)
	{
		has_decayed = true;
	}
	return has_decayed;
}

	/* Calculates if a position is eclipsed.  */
bool is_eclipsed(const double pos[3], const double sol[3], double *depth)
{
	double Rho[3], earth[3];

	/* Determine partial eclipse */
	double sd_earth = ArcSin(xkmper / vec3_length(pos));
	vec3_sub(sol, pos, Rho);
	double sd_sun = ArcSin(sr / vec3_length(Rho));
	vec3_mul_scalar(pos, -1, earth);
	
	double delta = ArcCos( vec3_dot(sol, earth) / vec3_length(sol) / vec3_length(earth) );
	*depth = sd_earth - sd_sun - delta;

	if (sd_earth < sd_sun) return false;
	else if (*depth >= 0) return true;
	else return false;
}

bool predict_is_eclipsed(const predict_orbit_t *orbit)
{
	return orbit->eclipsed;
}

double predict_eclipse_depth(const predict_orbit_t *orbit)
{
	return orbit->eclipse_depth;
}

double predict_squint_angle(const predict_observer_t *observer, const predict_orbit_t *orbit, double alon, double alat)
{
	double squint;
	if (orbit->ephemeris == EPHEMERIS_SDP4) {
		const struct _sdp4* sdp4 = (struct _sdp4*)orbit->ephemeris_data;

		double bx = cos(alat)*cos(alon + sdp4->deep_arg.omgadf);
		double by = cos(alat)*sin(alon + sdp4->deep_arg.omgadf);
		double bz = sin(alat);

		double cx = bx;
		double cy = by*cos(sdp4->xinck) - bz*sin(sdp4->xinck);
		double cz = by*sin(sdp4->xinck) + bz*cos(sdp4->xinck);
		double ax = cx*cos(sdp4->xnodek) - cy*sin(sdp4->xnodek);
		double ay = cx*sin(sdp4->xnodek) + cy*cos(sdp4->xnodek);
		double az = cz;

		struct predict_observation obs;
		predict_observe_orbit(observer, orbit, &obs);
		squint = acos(-(ax*obs.range_x + ay*obs.range_y + az*obs.range_z)/obs.range);

	} else {
		squint = nan("");
	}
	return squint;
}
