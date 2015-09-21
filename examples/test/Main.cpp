#include <stdio.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <limits>
#include <predict/predict.h>

#include "TestCase.h"

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


bool fuzzyCompare(const double &x, const double &y, const double &epsilon = std::numeric_limits<double>::epsilon())
{
	return fabs(x - y) < epsilon;
}

int runtest(const char *filename)
{
	
	// Load testcase
	TestCase testcase;
	if (testcase.loadFromFile(filename)) {
		fprintf(stderr, "Failed to load testfile: %s\n", filename);
		return -1;
	}

	// Get TLE
	char *tle[2];
	testcase.getTLE(tle);

	// Create orbit object
	predict_orbit_t *orbit = predict_create_orbit((const char**)tle);
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

		// Predict
		predict_orbit(orbit, predict_to_julian(time));

		// Observe
		struct predict_observation orbit_obs;
		predict_observe_orbit(obs, orbit, &orbit_obs);
		
		// Compare values
		string failed = "";
		if (!fuzzyCompare(orbit->latitude*180.0/M_PI, lat, 1E-8)) {
			failed += "(latitude)";
		}
		if (!fuzzyCompare(orbit->longitude*180.0/M_PI, lon, 1E-8)) {
			failed += "(longitude)";
		}
		if (!fuzzyCompare(orbit->altitude, alt, 1E-3)) {
			failed += "(altitude)";
		}
		if (!fuzzyCompare(orbit_obs.azimuth*180.0/M_PI, az, 1E-2)) {
			failed += "(azimuth)";
		}
		if (!fuzzyCompare(orbit_obs.elevation*180.0/M_PI, el, 1E-2)) {
			failed += "(elevation)";
		}

		// Failed?
		if (failed != "") {
			
			cout << filename << ": failed at data line " << line << ": " << failed << endl;

			printf("%.8f, %.8f/%.8f, %.8f/%.8f, %.3f/%.3f, %.3f/%.3f, %.3f/%.3f\n", time,
					orbit->latitude*180.0/M_PI, lat, 
					orbit->longitude*180.0/M_PI, lon,
					orbit->altitude, alt,
					orbit_obs.azimuth*180.0/M_PI, az,
					orbit_obs.elevation*180.0/M_PI, el);

			retval = -1;
		}

		// Increment data line number
		++line;

	}

	return retval;

}
