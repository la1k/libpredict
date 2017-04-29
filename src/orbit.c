#define _XOPEN_SOURCE 600
#include <math.h>
#include <string.h>
#include "defs.h"
#include "unsorted.h"
#include "sdp4.h"
#include "sgp4.h"
#include "sun.h"

bool is_eclipsed(const double pos[3], const double sol[3], double *depth);
bool predict_decayed(const predict_orbital_elements_t *orbital_elements, predict_julian_date_t time);

predict_orbital_elements_t* predict_parse_tle(const char *tle_line_1, const char *tle_line_2)
{
	double tempnum;
	predict_orbital_elements_t *m = (predict_orbital_elements_t*)malloc(sizeof(predict_orbital_elements_t));
	if (m == NULL) return NULL;

	m->satellite_number = atol(SubString(tle_line_1,2,6));
	m->element_number = atol(SubString(tle_line_1,64,67));
	m->epoch_year = atoi(SubString(tle_line_1,18,19));
	strncpy(m->designator, SubString(tle_line_1,9,16),8);
	m->epoch_day = atof(SubString(tle_line_1,20,31));
	m->inclination = atof(SubString(tle_line_2,8,15));
	m->right_ascension = atof(SubString(tle_line_2,17,24));
	m->eccentricity = 1.0e-07*atof(SubString(tle_line_2,26,32));
	m->argument_of_perigee = atof(SubString(tle_line_2,34,41));
	m->mean_anomaly = atof(SubString(tle_line_2,43,50));
	m->mean_motion = atof(SubString(tle_line_2,52,62));
	m->derivative_mean_motion  = atof(SubString(tle_line_1,33,42));
	tempnum=1.0e-5*atof(SubString(tle_line_1,44,49));
	m->second_derivative_mean_motion = tempnum/pow(10.0,(tle_line_1[51]-'0'));
	tempnum=1.0e-5*atof(SubString(tle_line_1,53,58));
	m->bstar_drag_term = tempnum/pow(10.0,(tle_line_1[60]-'0'));
	m->revolutions_at_epoch = atof(SubString(tle_line_2,63,67));

	/* Period > 225 minutes is deep space */
	double ao, xnodp, dd1, dd2, delo, a1, del1, r1;
	double temp = twopi/xmnpda/xmnpda;
	double xno = m->mean_motion*temp*xmnpda; //from old TLE struct
	dd1=(xke/xno);
	dd2=tothrd;
	a1=pow(dd1,dd2);
	r1=cos(m->inclination*M_PI/180.0);
	dd1=(1.0-m->eccentricity*m->eccentricity);
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
			predict_destroy_orbital_elements(m);
			return NULL;
		}
		// Initialize ephemeris data structure
		sdp4_init(m, (struct _sdp4*)m->ephemeris_data);

	} else {
		m->ephemeris = EPHEMERIS_SGP4;
		
		// Allocate memory for ephemeris data
		m->ephemeris_data = malloc(sizeof(struct _sgp4));

		if (m->ephemeris_data == NULL) {
			predict_destroy_orbital_elements(m);
			return NULL;
		}
		// Initialize ephemeris data structure
		sgp4_init(m, (struct _sgp4*)m->ephemeris_data);
	}

	return m;
}

void predict_destroy_orbital_elements(predict_orbital_elements_t *m)
{
	if (m == NULL) return;

	if (m->ephemeris_data != NULL) {
		free(m->ephemeris_data);
	}

	free(m);
}

bool predict_is_geosynchronous(const predict_orbital_elements_t *m)
{
	return (m->mean_motion >= GEOSYNCHRONOUS_LOWER_MEAN_MOTION)
		&& (m->mean_motion <= GEOSYNCHRONOUS_UPPER_MEAN_MOTION)
		&& (fabs(m->eccentricity) <= GEOSYNCHRONOUS_ECCENTRICITY_THRESHOLD)
		&& (fabs(m->inclination) <= GEOSYNCHRONOUS_INCLINATION_THRESHOLD_DEGREES);
}

double predict_apogee(const predict_orbital_elements_t *m)
{
	double sma = 331.25*exp(log(1440.0/m->mean_motion)*(2.0/3.0));
	return sma*(1.0+m->eccentricity)-xkmper;
}
		
double predict_perigee(const predict_orbital_elements_t *m)
{
	double xno = m->mean_motion*twopi/xmnpda;
	double a1=pow(xke/xno,tothrd);
	double cosio=cos(m->inclination*M_PI/180.0);
	double theta2=cosio*cosio;
	double x3thm1=3*theta2-1.0;
	double eosq=m->eccentricity*m->eccentricity;
	double betao2=1.0-eosq;
	double betao=sqrt(betao2);
	double del1=1.5*ck2*x3thm1/(a1*a1*betao*betao2);
	double ao=a1*(1.0-del1*(0.5*tothrd+del1*(1.0+134.0/81.0*del1)));
	double delo=1.5*ck2*x3thm1/(ao*ao*betao*betao2);
	double aodp=ao/(1.0-delo);

	return (aodp*(1-m->eccentricity)-ae)*xkmper;
}

bool predict_aos_happens(const predict_orbital_elements_t *m, double latitude)
{
	/* This function returns true if the satellite pointed to by
	   "x" can ever rise above the horizon of the ground station. */

	double lin, apogee;

	if (m->mean_motion==0.0)
		return false;
	else
	{
		lin = m->inclination;

		if (lin >= 90.0) lin = 180.0-lin;

		apogee = predict_apogee(m);

		if ((acos(xkmper/(apogee+xkmper))+(lin*M_PI/180.0)) > fabs(latitude))
			return true;
		else
			return false;
	}
}

/* This is the stuff we need to do repetitively while tracking. */
/* This is the old Calc() function. */
int predict_orbit(const predict_orbital_elements_t *orbital_elements, struct predict_orbit *m, double utc)
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
	double epoch = 1000.0*orbital_elements->epoch_year + orbital_elements->epoch_day;
	double jul_epoch = Julian_Date_of_Epoch(epoch);
	double tsince = (julTime - jul_epoch)*xmnpda;

	/* Call NORAD routines according to deep-space flag. */
	struct model_output output;
	switch (orbital_elements->ephemeris) {
		case EPHEMERIS_SDP4:
			sdp4_predict((struct _sdp4*)orbital_elements->ephemeris_data, tsince, &output);
			break;
		case EPHEMERIS_SGP4:
			sgp4_predict((struct _sgp4*)orbital_elements->ephemeris_data, tsince, &output);
			break;
		default:
			//Panic!
			return -1;
	}
	m->position[0] = output.pos[0];
	m->position[1] = output.pos[1];
	m->position[2] = output.pos[2];
	m->velocity[0] = output.vel[0];
	m->velocity[1] = output.vel[1];
	m->velocity[2] = output.vel[2];
	m->phase = output.phase;
	m->argument_of_perigee = output.omgadf;
	m->inclination = output.xinck;
	m->right_ascension = output.xnodek;

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
	double xno = orbital_elements->mean_motion*temp*xmnpda;
	double xmo = orbital_elements->mean_anomaly * M_PI / 180.0;
	m->revolutions = (long)floor((xno*xmnpda/(M_PI*2.0) + age*orbital_elements->bstar_drag_term)*age + xmo/(2.0*M_PI)) + orbital_elements->revolutions_at_epoch;

	//calculate whether orbit is decayed
	m->decayed = predict_decayed(orbital_elements, utc);

	return 0;
}

bool predict_decayed(const predict_orbital_elements_t *orbital_elements, predict_julian_date_t time)
{
	double satepoch;
	satepoch=DayNum(1,0,orbital_elements->epoch_year)+orbital_elements->epoch_day;

	bool has_decayed = false;
	if (satepoch + ((16.666666 - orbital_elements->mean_motion)/(10.0*fabs(orbital_elements->derivative_mean_motion))) < time)
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

bool predict_is_eclipsed(const struct predict_orbit *orbit)
{
	return orbit->eclipsed;
}

double predict_eclipse_depth(const struct predict_orbit *orbit)
{
	return orbit->eclipse_depth;
}

double predict_squint_angle(const predict_observer_t *observer, const struct predict_orbit *orbit, double alon, double alat)
{
	double bx = cos(alat)*cos(alon + orbit->argument_of_perigee);
	double by = cos(alat)*sin(alon + orbit->argument_of_perigee);
	double bz = sin(alat);

	double cx = bx;
	double cy = by*cos(orbit->inclination) - bz*sin(orbit->inclination);
	double cz = by*sin(orbit->inclination) + bz*cos(orbit->inclination);
	double ax = cx*cos(orbit->right_ascension) - cy*sin(orbit->right_ascension);
	double ay = cx*sin(orbit->right_ascension) + cy*cos(orbit->right_ascension);
	double az = cz;

	struct predict_observation obs;
	predict_observe_orbit(observer, orbit, &obs);
	double squint = acos(-(ax*obs.range_x + ay*obs.range_y + az*obs.range_z)/obs.range);

	return squint;
}
