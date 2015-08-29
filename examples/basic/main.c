#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include <predict/predict.h>

int main(int argc, char **argv)
{

	const char *tle[2] = {
		"1 25544U 98067A   15129.86961041  .00015753  00000-0  23097-3 0  9998",
		"2 25544  51.6464 275.3867 0006524 289.1638 208.5861 15.55704207942078"};

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
		predict_julian_date_t curr_time = predict_to_julian(time(NULL));
	
		// Predict ISS
		orbit_predict(iss, curr_time);
		printf("ISS: lat=%f, lon=%f, alt=%f, eclipsed=%i (%.2f)\n", iss->latitude*180.0/M_PI, iss->longitude*180.0/M_PI, 
				iss->altitude, orbit_is_eclipsed(iss), orbit_eclipse_depth(iss)*180.0/M_PI);
	
		// Observe ISS
		struct observation iss_obs;
		observer_find_orbit(obs, iss, &iss_obs);
		printf("ISS: %f (rate: %f), %f (rate: %f)\n", iss_obs.azimuth*180.0/M_PI, iss_obs.azimuth_rate*180.0/M_PI, iss_obs.elevation*180.0/M_PI, iss_obs.elevation_rate*180.0/M_PI);

		// Apparent elevation
		double apparent_elevation = predict_apparent_elevation(iss_obs.elevation);
		printf("Apparent ISS elevation: %.2f\n", apparent_elevation*180.0/M_PI);

		// Predict and observe MOON
		struct observation moon_obs;
		observer_find_moon(obs, curr_time, &moon_obs);
		printf("MOON: %f, %f\n", moon_obs.azimuth*180.0/M_PI, moon_obs.elevation*180.0/M_PI);

		// Predict and observe SUN
		struct observation sun_obs;
		observer_find_sun(obs, curr_time, &sun_obs);
		printf("SUN: %f, %f\n", sun_obs.azimuth*180.0/M_PI, sun_obs.elevation*180.0/M_PI);

		//Sleep
		usleep(1000000);
	}

	// Free memory
	orbit_destroy(iss);
	observer_destroy(obs);

	return 0;
}
