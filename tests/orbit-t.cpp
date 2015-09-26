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

	// Create orbit objects

	// Used in lower bound in value check
	predict_orbit_t *orbit_lower = predict_create_orbit(tle);
	if (!orbit_lower) {
		fprintf(stderr, "Failed to initialize orbit from tle!");
		return -1;
	}

	// Used in upper bound in value check
	predict_orbit_t *orbit_upper = predict_create_orbit(tle);
	if (!orbit_upper) {
		fprintf(stderr, "Failed to initialize orbit from tle!");
		return -1;
	}

	// Create observer object
	predict_observer_t *obs = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		return -1;
	}

	bool check_squint_angle = testcase.containsValidAlonAlat() && (orbit_lower->ephemeris == EPHEMERIS_SDP4);
	
	// Test
	int retval = 0;
	int line = 1;
	for (auto d : testcase.data()) {
		double time = d[0];
		double lat = d[1];
		double lon = d[2];
		double alt = d[3];
		double az = d[4];
		double el = d[5];
		double squint = d[7];

		// Compare values within (time - 1, time + 1) (i.e. up time + 1, but not including time + 1)
		// (since we don't know the exact time predict generated its data, only within an error of 1 second)
		const int DIFF = 1;
		struct predict_observation orbit_obs_lower;
		struct predict_observation orbit_obs_upper;

		// Lower bound
		predict_orbit(orbit_lower, predict_to_julian(time));
		predict_observe_orbit(obs, orbit_lower, &orbit_obs_lower);

		// Upper bound
		predict_orbit(orbit_upper, predict_to_julian(time + DIFF));
		predict_observe_orbit(obs, orbit_upper, &orbit_obs_upper);

		// Check values
		string failed = "";
		if (!fuzzyCompareWithBoundaries(orbit_lower->latitude*180.0/M_PI, orbit_upper->latitude*180/M_PI, lat)) {
			failed += "(latitude)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_lower->longitude*180.0/M_PI, orbit_upper->longitude*180/M_PI, lon)) {
			failed += "(longitude)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_lower->altitude, orbit_upper->altitude, alt)) {
			failed += "(altitude)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_obs_lower.azimuth*180.0/M_PI, orbit_obs_upper.azimuth*180.0/M_PI, az)) {
			failed += "(azimuth)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_obs_lower.elevation*180.0/M_PI, orbit_obs_upper.elevation*180.0/M_PI, el)) {
			failed += "(elevation)";
		}

		double squint_angle_lower, squint_angle_upper;
		if (check_squint_angle) {
			squint_angle_lower = predict_squint_angle(obs, orbit_lower, testcase.alon()*M_PI/180.0, testcase.alat()*M_PI/180.0)*180.0/M_PI;
			squint_angle_upper = predict_squint_angle(obs, orbit_upper, testcase.alon()*M_PI/180.0, testcase.alat()*M_PI/180.0)*180/M_PI;
			if (!fuzzyCompareWithBoundaries(squint_angle_lower, squint_angle_upper, squint)) {
				failed += "(squint)";
			}
		}

		// Failed?
		if (failed != "") {
			cout << filename << ": failed at data line " << line << ": " << failed << endl;

			printf("%.8f, %.8f/%.8f/%.8f, %.8f/%.8f/%.8f, %.3f/%.3f/%.3f, %.3f/%.3f/%.3f, %.3f/%.3f/%.3f", time,
					orbit_lower->latitude*180.0/M_PI, lat, orbit_upper->latitude*180.0/M_PI,
					orbit_lower->longitude*180.0/M_PI, lon, orbit_upper->longitude*180.0/M_PI,
					orbit_lower->altitude, alt, orbit_upper->altitude,
					orbit_obs_lower.azimuth*180.0/M_PI, az, orbit_obs_upper.azimuth*180.0/M_PI,
					orbit_obs_lower.elevation*180.0/M_PI, el, orbit_obs_upper.elevation*180.0/M_PI);
			if (check_squint_angle) {
				printf(" %.8f/%.8f/%.8f", squint_angle_lower, squint, squint_angle_upper);
			}
			printf("\n");

			retval = -1;
		}

		// Increment data line number
		++line;
	}

	return retval;
}
