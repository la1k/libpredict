#ifndef ORBIT_H_
#define ORBIT_H_

#include <predict/tle.h>
#include <stdbool.h>

typedef enum {
  EPHEMERIS_SGP4 = 0,
  EPHEMERIS_SDP4 = 1,
  EPHEMERIS_SGP8 = 2,
  EPHEMERIS_SDP8 = 3
} ephemeris_t;

struct orbit {
	char name[128];
	double time;
	double position[3];
	double velocity[3];

	double latitude;
	double longitude;
	double altitude;
	int eclipsed;
	ephemeris_t ephemeris;
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

	///Ephemeris data structure pointer
	void *ephemeris_data;

};

typedef struct orbit orbit_t;

orbit_t *orbit_create(const char *tle[]);
void orbit_destroy(orbit_t *orbit);

int orbit_predict(orbit_t *x, double utc);
bool orbit_is_geostationary(const orbit_t *x);
double orbit_apogee(const orbit_t *x);
double orbit_perigee(const orbit_t *x);
bool orbit_aos_happens(const orbit_t *x, double latitude);

/* return true if orbit has decayed */
bool orbit_decayed(const orbit_t *x);


#endif
