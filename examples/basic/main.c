#include <stdio.h>
#include "orbit.h"

int main(int argc, char **argv)
{

	const char *tle[2] = {
		"1 25544U 98067A   12317.22408565  .00013449  00000-0  23146-3 0  6331", 
		"2 25544  51.6464 106.0986 0015548 274.8897  81.8006 15.51295927801004" };

	struct orbit iss;
	orbit_init(&iss, tle);

	for (;;) {
		
		//Get current time:
		double time = CurrentDaynum();

		orbit_predict(&iss, time);

		printf("%f, %f, %f\n", iss.latitude*180.0/M_PI, iss.longitude*180.0/M_PI, iss.altitude);
		
		//Sleep
		usleep(100000);

	}

	return 0;
}
