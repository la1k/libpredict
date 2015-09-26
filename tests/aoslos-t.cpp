#include <stdio.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <predict/predict.h>
#include "testcase_reader.h"
#include <iostream>
using namespace std;

int runtest(const char *filename);

int main(int argc, char **argv)
{
	// Check arguments
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " <testfiles>" << endl;
		return -1;
	}

	// Test all provided test files
	int retval = 0;
	for (int i = 1; i < argc; ++i) {
		if (runtest(argv[i])) {
			cout << argv[i] << ": failed" << endl;
			retval = -1;
		} else {
			cout << argv[i] << ": OK" << endl;
		}
	}

	return retval;
}

int runtest(const char *filename)
{
	// Load testcase
	TestCaseReader testcase;
	testcase.loadFromFile(filename);
	if (!(testcase.containsValidData() && (testcase.containsValidQth()) && (testcase.containsValidTLE()))) {
		fprintf(stderr, "Failed to load testfile: %s\n", filename);
		return -1;
	}

	// Get TLE
	char *tle[2];
	testcase.getTLE(tle);

	// Create orbit object
	predict_orbit_t *orbit = predict_create_orbit(tle);
	if (!orbit) {
		fprintf(stderr, "Failed to initialize orbit from tle!");
		return -1;
	}

	// Create observer object
	predict_observer_t *obs = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		return -1;
	}
	
	int retval = 0;
	int line = 1;

	// Use first available time as start time for AOS/LOS finding
	double start_time = testcase.data()[0][0];

	predict_julian_date_t next_aos_time = predict_next_aos(obs, orbit, predict_to_julian(start_time));
	predict_julian_date_t next_los_time = predict_next_los(obs, orbit, predict_to_julian(start_time)); //can be LOS of current pass, if satellite is in range
	
	double time_diff = 1.0/(60.0*60.0*24.0); //1 second
	predict_julian_date_t curr_time = predict_to_julian(start_time);

	const double ELEVATION_THRESH = 0.3; //this elevation threshold was used internal AOS/LOS calculations, so should hold also here.

	// Check times until the AOS
	while (curr_time < next_aos_time) {
		struct predict_observation orbit_obs;
		predict_orbit(orbit, curr_time);
		predict_observe_orbit(obs, orbit, &orbit_obs);

		if ((next_los_time < next_aos_time) && (curr_time < next_los_time)) {
			//satellite should be above the horizon (satellite was in range at the start time)
			if (orbit_obs.elevation*180/M_PI < -ELEVATION_THRESH) {
				fprintf(stderr, "AOS failed sanity check: Satellite currently in range according to LOS/AOS state, but not above the horizon. Elevation: %f\n", orbit_obs.elevation*180/M_PI);
				return -1;
			}
		} else if (curr_time < next_aos_time) {
			//satellite should be below the horizon
			if (orbit_obs.elevation*180/M_PI > ELEVATION_THRESH){
				fprintf(stderr, "AOS failed sanity check: Satellite has not reached AOS, but is above the horizon. Elevation: %f\n", orbit_obs.elevation*180/M_PI);
				return -1;
			}
		} 
		curr_time += time_diff;
	}

	// Check times up till LOS
	next_los_time = predict_next_los(obs, orbit, predict_to_julian(curr_time)); //recalculating within the pass in case LOS was for previous pass
	while (curr_time < next_los_time) {
		struct predict_observation orbit_obs;
		predict_orbit(orbit, curr_time);
		predict_observe_orbit(obs, orbit, &orbit_obs);

		//satellite should be above the horizon
		if (orbit_obs.elevation*180/M_PI < -ELEVATION_THRESH) {
			fprintf(stderr, "AOS/LOS failed sanity check: Satellite has reached AOS, waiting for LOS, but not above the horizon. Elevation: %f\n", orbit_obs.elevation*180/M_PI);
			return -1;
		}
		
		
		curr_time += time_diff;
	}

	return retval;
}
