#ifndef _OBSERVER_H_
#define _OBSERVER_H_

#include <predict/orbit.h>

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

struct observation {
	///UTC time                
	double time;                       
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

observer_t *observer_create(const char *name, double lat, double lon, double alt);
void observer_destroy(observer_t *obs);

void observer_find_orbit(const observer_t *observer, const orbit_t *orbit, struct observation *obs);
void observer_find_moon(const observer_t *observer, double time, struct observation *obs);
void observer_find_sun(const observer_t *observer, double time, struct observation *obs);

/* Find next AOS from time start_utc (ignore previous AOS of current pass if satellite is in range) */
double observer_get_next_aos(const observer_t *observer, orbit_t *orbit, double start_utc);

/* Find next LOS from time start_utc (LOS of an upcoming pass or the current pass if satellite is in range) */
double observer_get_next_los(const observer_t *observer, orbit_t *orbit, double start_utc);

/* calculate doppler shift of a given downlink frequency with respect to the observer */
double observer_get_doppler_shift(const observer_t *observer, const orbit_t *orbit, double frequency);

#endif // ifndef _OBSERVER_H_
