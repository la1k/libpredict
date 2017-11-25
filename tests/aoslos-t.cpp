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

/**
 * Predict AOS/LOS for several passes after the start time, and check that the times are consistent with each other, that elevation is non-zero in the center
 * and that the elevation rates are correct at each point.
 *
 * \param orbital_elements Orbital elements
 * \param observer Observer
 * \param start_time Start time
 * \return 0 on success, -1 on failure
 **/
int aoslos_timepoint_consistency_test(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, double start_time);

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
	predict_orbital_elements_t *orbital_elements = predict_parse_tle(tle[0], tle[1]);

	// Create observer object
	predict_observer_t *obs = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		return -1;
	}

	// Use first available time as start time for AOS/LOS finding
	double start_time = predict_to_julian(testcase.data()[0][0]);

	//check whether the pass makes sense wrt elevation and/or elevation rate at start, end and middle of pass
	if (aoslos_timepoint_consistency_test(orbital_elements, obs, start_time) != 0) {
		return -1;
	}

	return 0;
}

//The tolerance used for fine-tuning the elevations in the aoslos-functions
#define ELEVATION_ZERO_TOLERANCE 0.3

//Expected time between two passes should at least be 20 minutes
#define PASS_TIME_THRESH (20.0/(24*60))

int aoslos_timepoint_consistency_test(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, double start_time)
{
	if (predict_is_geosynchronous(orbital_elements)) {
		return 0;
	}

	struct predict_observation aos, los;
	aos = predict_next_aos(observer, orbital_elements, start_time);
	double aos_time = aos.time;
	los = predict_next_los(observer, orbital_elements, aos_time);
	double los_time = los.time;

	struct predict_position orbit;
	struct predict_observation observation;

	const int NUM_PASSES = 10;
	for (int i=0; i < NUM_PASSES; i++) {
		if (los_time <= aos_time) {
			fprintf(stderr, "los time not strictly larger than aos time: %f %f\n", aos_time, los_time);
			return -1;
		}

		predict_orbit(orbital_elements, &orbit, aos_time);
		predict_observe_orbit(observer, &orbit, &observation);
		double elevation_rate_aos = observation.elevation_rate;

		predict_orbit(orbital_elements, &orbit, los_time);
		predict_observe_orbit(observer, &orbit, &observation);
		double elevation_rate_los = observation.elevation_rate;

		if ((elevation_rate_aos <= 0) || (elevation_rate_los >= 0)) {
			fprintf(stderr, "Elevation rates do not have correct signs for aos and los times: %f (%f) %f (%f)\n", elevation_rate_aos, aos_time, elevation_rate_los, los_time);
			return -1;
		}

		double midpass = (los_time - aos_time)/2.0 + aos_time;
		predict_orbit(orbital_elements, &orbit, midpass);
		predict_observe_orbit(observer, &orbit, &observation);

		if (observation.elevation <= 0.0) {
			fprintf(stderr, "Elevation is negative during the middle of a pass: %f %f %f\n", aos_time, midpass, los_time);
			return -1;
		}

		double prev_los_time = los_time;
		double prev_aos_time = aos_time;
		struct predict_observation aos, los;
		aos = predict_next_aos(observer, orbital_elements, los_time);
		aos_time = aos.time;
		los = predict_next_los(observer, orbital_elements, aos_time);
		los_time = los.time;

		if (!((los_time > aos_time) && (los_time > prev_los_time) && (los_time > prev_aos_time) && (aos_time > prev_los_time))) {
			fprintf(stderr, "New AOS/LOS not strictly larger than previous AOS/LOS\n");
			return -1;
		}

		if ((aos_time - prev_los_time) < PASS_TIME_THRESH) {
			fprintf(stderr, "Time between passes not significantly different: %f minutes\n", (aos_time - prev_los_time)*24*60);
			return -1;
		}

		//check that time between LOS time and next AOS time produces negative elevations
		double timestep = 1.0/(24.0*50);
		double curr_time = prev_los_time + timestep;
		while (curr_time < aos_time - timestep) {
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &observation);
			if (observation.elevation*180.0/M_PI >= ELEVATION_ZERO_TOLERANCE) {
				fprintf(stderr, "AOS/LOS failed consistency check between two passes: %f %f %f %f\n", curr_time, observation.elevation, prev_los_time, aos_time);
				return -1;
			}
			curr_time += timestep;
		}
	}

	return 0;
}
