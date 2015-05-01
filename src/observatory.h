#ifndef _OBSERVATORY_H_
#define _OBSERVATORY_H_

struct observatory {
	///Observatory name
	char name[128];
	///Latitude (WGS84, radians)
	double latitude;
	///Longitude (WGS84, radians)
	double longitude;
	///Altitude (WGS84, meters)
	double altitude;
};

#endif // ifndef _OBSERVATORY_H_
