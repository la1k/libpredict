#ifndef ORBIT_H_
#define ORBIT_H_

#include "vec3.h"
#include "sdp4.h"
#include "sgp4.h"

#define EPHEMERIS_SGP4 		0
#define EPHEMERIS_SDP4		1
#define EPHEMERIS_SGP8		2
#define EPHEMERIS_SDP8		3

struct orbit {
	char name[128];
	double time;
	double position[3];
	double velocity[3];

	double latitude;
	double longitude;
	double altitude;
	int eclipsed;
	int ephemeris;
	///Original TLE line number one:
	char line1[70];
	///Original TLE line number two:
	char line2[70];
	///Original tle_t used to hold processed tle parameters used in calculations.
	tle_t tle;

	///Satellite number (line 1, field 2)
	long catnum;
	///Element number (line 1, field 13)
	long setnum;
	///International designator (line 1, fields 4, 5, 6)
	char designator[10];
	///Epoch year (last two digits) (line 1, field 7)
	int year;
	///Epoch day (day of year and fractional portion of day, line 1, field 8)
	double refepoch;
	///Inclination (line 2, field 3)
	double incl;
	///Right Ascension of the Ascending Node [Degrees] (line 2, field 4)
	double raan;
	///Eccentricity (decimal point assumed) (line 2, field 5)
	double eccn;
	///Argument of Perigee [Degrees] (line 2, field 6)
	double argper;
	///Mean Anomaly [Degrees] (line 2, field 7)
	double meanan;
	///Mean Motion [Revs per day] (line 2, field 8)
	double meanmo;
	///First Time Derivative of the Mean Motion divided by two (line 1, field 9)
	double drag;
	///Second Time Derivative of Mean Motion divided by six (decimal point assumed, line 1, field 10)
	double nddot6;
	///BSTAR drag term (decimal point assumed, line 1, field 11)
	double bstar;
	///Orbital number (line 2, field 9)
	long orbitnum;

	struct _sdp4 sdp4;
	struct _sgp4 sgp4;

};

int orbit_init(struct orbit *x, const char *tle[]);
int orbit_predict(struct orbit *x, double utc);
int orbit_is_geostationary(const struct orbit *x);
double orbit_apogee(const struct orbit *x);
double orbit_perigee(const struct orbit *x);
int orbit_aos_happens(const struct orbit *x, double latitude);

#endif
