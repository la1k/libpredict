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
void observer_find_orbit(const observer_t *observer, const orbit_t *orbit, struct observation *obs);
void observer_find_moon(const observer_t *observer, double time, struct observation *obs);
void observer_find_sun(const observer_t *observer, double time);

#endif // ifndef _OBSERVER_H_
