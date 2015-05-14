#ifndef _PREDICT_H_
#define _PREDICT_H_

#include <time.h>
#include <stdbool.h>

/**
 * The representation of time used by libpredict: The number of days since 31Dec79 00:00:00 UTC. 
 **/
typedef double predict_julian_date_t;

/**
 * Convert time_t in UTC to Julian date in UTC.
 *
 * \param time Time in UTC
 * \return Julian day in UTC
 **/
predict_julian_date_t predict_get_julian_date_from_time(time_t time);

/**
 * Convert Julian date in UTC back to a time_t in UTC. 
 *
 * \param date Julian date in UTC
 * \return Time in UTC
 **/
time_t predict_get_time_from_julian_date(predict_julian_date_t date);

/**
 * Container for processed TLE data from TLE strings.
 **/
typedef struct	{
	double epoch;
	double xndt2o;
	double xndd6o;
	double bstar;
	double xincl;
	double xnodeo;
	double eo;
	double omegao;
	double xmo;
	double xno;
 	int catnr;
	int elset;
	int revnum;
}tle_t; 

/**
 * Simplified perturbation models used in modeling the satellite orbits. 
 **/
typedef enum {
  EPHEMERIS_SGP4 = 0,
  EPHEMERIS_SDP4 = 1,
  EPHEMERIS_SGP8 = 2,
  EPHEMERIS_SDP8 = 3
} ephemeris_t;

/**
 * Satellite orbit definitions, according to defined NORAD TLE. 
 **/
struct orbit {
	///Name of satellite
	char name[128];

	///Timestamp for last call to orbit_predict
	predict_julian_date_t time;
	///ECI position in km
	double position[3];
	///ECI velocity in km/s
	double velocity[3];

	///Latitude in radians, northing/easting
	double latitude;
	///Longitude in radians, northing/easting
	double longitude;
	///Altitude in meters
	double altitude;
	///Whether satellite is eclipsed by the earth
	int eclipsed;
	///Which perturbation model to use
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

/**
 * Create orbit structure. 
 * \param tle NORAD two-line element set as string array
 * \return Allocated orbit structure
 **/
orbit_t *orbit_create(const char *tle[]);

/**
 * Free memory allocated in orbit structure. 
 * \param orbit Orbit to free
 **/
void orbit_destroy(orbit_t *orbit);

/**
 * Main prediction function. Predict satellite orbit at given time. 
 * \param x Satellite orbit
 * \param time Julian day in UTC
 * \return 0 if everything went fine
 **/
int orbit_predict(orbit_t *x, predict_julian_date_t time);

/**
 * Find whether orbit is geostationary. 
 *
 * \param x Satellite orbit
 * \return true if orbit is geostationary, otherwise false
 **/
bool orbit_is_geostationary(const orbit_t *x);

/** 
 * Get apogee of satellite orbit. 
 *
 * \param x Satellite orbit
 * \return Apogee of orbit
 **/
double orbit_apogee(const orbit_t *x);

/**
 * Get perigee of satellite orbit. 
 *
 * \param x Satellite orbit
 * \return Perigee of orbit
 **/
double orbit_perigee(const orbit_t *x);

/**
 * Find whether an AOS can ever happen on the given latitude. 
 *
 * \param x Satellite orbit
 * \param latitude Latitude of ground station
 * \return true if AOS can happen, otherwise false
 **/
bool orbit_aos_happens(const orbit_t *x, double latitude);

/** 
 * Find whether an orbit has decayed.
 *
 * \param x Current state of orbit
 * \return true if orbit has decayed, otherwise false
 **/
bool orbit_decayed(const orbit_t *x);

/**
 * Observation point/ground station (QTH).
 **/
typedef struct observer {
	///Observatory name
	char name[128];
	///Latitude (WGS84, radians)
	double latitude;
	///Longitude (WGS84, radians)
	double longitude;
	///Altitude (WGS84, meters)
	double altitude;
} observer_t;

/**
 * Data relevant for a relative observation of an orbit or similar with respect to an observation point.
 **/
struct observation {
	///UTC time                
	predict_julian_date_t time;                       
	///Azimut angle      
	double azimut; 
	///Elevation angle                               
	double elevation; 
	///Corrected elevation    
	double correctedElevation;         
	///Range (km) 
	double range;                        
	///Range vector                    
	double xRange, yRange, zRange; 
	///Range velocity (km/s) 
	double rangeDot;                    
	///Visible? 
	int visible; 
};

/**
 * Create observation point (QTH).
 *
 * \param name Name of observation point
 * \param lat Latitude in radians (easting/northing)
 * \param lon Longitude in radians (easting/northing)
 * \param alt Altitude in meters
 * \return Allocated observation point
 **/
observer_t *observer_create(const char *name, double lat, double lon, double alt);

/** 
 * Free observer.
 *
 * \param obs Observer to be freed.
 **/
void observer_destroy(observer_t *obs);

/** 
 * Find relative position of satellite with respect to an observer.
 *
 * \param observer Point of observation
 * \param orbit Satellite orbit
 * \param obs Return of object for position of the satellite relative to the observer.
 **/
void observer_find_orbit(const observer_t *observer, const orbit_t *orbit, struct observation *obs);

/**
 * Estimate relative position of the moon.
 *
 * \param observer Point of observation
 * \param time Time of observation
 * \param obs Return object for position of the moon relative to the observer
 **/
void observer_find_moon(const observer_t *observer, predict_julian_date_t time, struct observation *obs);

/** 
 * Estimate relative position of the sun.
 *
 * \param observer Point of observation
 * \param time Time of observation
 * \param obs Return object for position of the sun relative to the observer
 **/
void observer_find_sun(const observer_t *observer, predict_julian_date_t time, struct observation *obs);

/** 
 * Find next acquisition of signal (AOS) of satellite (when the satellite rises above the horizon). Ignores previous AOS of current pass if the satellite is in range at the start time. 
 *
 * \param observer Point of observation
 * \param orbit Satellite orbit
 * \param start_time Start time for AOS search
 * \return Time of AOS
 **/
predict_julian_date_t observer_get_next_aos(const observer_t *observer, orbit_t *orbit, predict_julian_date_t start_time);

/** 
 * Find next loss of signal (LOS) of satellite (when the satellite goes below the horizon). Finds LOS of the current pass if the satellite currently is in range, finds LOS of next pass if not.
 *
 * \param observer Point of observation
 * \param orbit Satellite orbit
 * \param start_time Start time for LOS search
 * \return Time of LOS
 **/
predict_julian_date_t observer_get_next_los(const observer_t *observer, orbit_t *orbit, predict_julian_date_t start_time);

/**
 * Calculate doppler shift of a given downlink frequency with respect to the observer. 
 *
 * \param observer Point of observation
 * \param orbit Current state of satellite orbit
 * \param downlink_frequency Downlink frequency of the satellite
 * \return The frequency difference from the original frequency
 **/
double observer_get_doppler_shift(const observer_t *observer, const orbit_t *orbit, double downlink_frequency);

#endif //_PREDICT_H_
