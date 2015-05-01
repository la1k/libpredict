#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <predict/predict.h>

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

	// Predict orbit
	int i;
	for (i=0;i<100;++i) {
		
		//Get current time:
		double time = CurrentDaynum();

		orbit_predict(iss, time);

		printf("%f, %f, %f\n", iss->latitude*180.0/M_PI, iss->longitude*180.0/M_PI, iss->altitude);
		
		//Sleep
		usleep(100000);
	}

	// Free memory
	orbit_destroy(iss);

	return 0;
}
