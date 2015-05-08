#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <predict/predict.h>
#include <predict/observer.h>

int main(int argc, char **argv)
{

	const char *tle[2] = {
		"1 25544U 98067A   12317.22408565  .00013449  00000-0  23146-3 0  6331", 
		"2 25544  51.6464 106.0986 0015548 274.8897  81.8006 15.51295927801004" };

	// Create orbit object
	orbit_t *iss = orbit_create(tle);
	if (!iss) {
		fprintf(stderr, "Failed to initialize orbit from tle!");
		exit(1);
	}

	// Create observer object
	observer_t *obs = observer_create("Me", 63.9*M_PI/180.0, 10.9*M_PI/180.0, 0);
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		exit(1);
	}


	// Predict orbit
	int i;
	for (i=0;i<100;++i) {
		
		//Get current time:
		double time = CurrentDaynum();

		// Predict ISS
		orbit_predict(iss, time);
		printf("ISS: lat=%f, lon=%f, alt=%f\n", iss->latitude*180.0/M_PI, iss->longitude*180.0/M_PI, iss->altitude);
	
		// Observe ISS
		struct observation iss_obs;
		observer_find_orbit(obs, iss, &iss_obs);
		printf("ISS: %f, %f\n", iss_obs.azimut*180.0/M_PI, iss_obs.elevation*180.0/M_PI);

		// Predict and observe MOON
		struct observation moon_obs;
		observer_find_moon(obs, time, &moon_obs);
		printf("MOON: %f, %f\n", moon_obs.azimut*180.0/M_PI, moon_obs.elevation*180.0/M_PI);

		// Predict and observe SUN
		struct observation sun_obs;
		observer_find_sun(obs, time, &sun_obs);
		printf("SUN: %f, %f\n", sun_obs.azimut*180.0/M_PI, sun_obs.elevation*180.0/M_PI);
	

		//Sleep
		usleep(100000);
	}

	// Free memory
	orbit_destroy(iss);

	return 0;
}
