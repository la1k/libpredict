#include <stdio.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <limits>
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
	if (!(testcase.containsValidData() && (testcase.containsValidQth()))) {
		fprintf(stderr, "Failed to load testfile: %s\n", filename);
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
	for (size_t i=0; i < testcase.data().size(); i++) {
		std::vector<double> d = testcase.data()[i];
		double time = d[0] - 2444238.5; // XXX: difference to other test cases input time is JD
		double az = d[1];
		double el = d[2];

		struct predict_observation jupiter_obs;

		predict_observe_jupiter(obs, time, &jupiter_obs);

		// Check values
		string failed = "";
		if (!fuzzyCompare(jupiter_obs.azimuth*180.0/M_PI, az, 0.06)) {
			failed += "(azimuth)";
		}
		if (!fuzzyCompare(jupiter_obs.elevation*180.0/M_PI, el, 0.06)) {
			failed += "(elevation)";
		}

		// Failed?
		if (failed != "") {

			cout << filename << ": failed at data line " << line << ": " << failed << endl;

			printf("%.8f, az_ref %.8f az %.8f diff %.8f, el_ref %.8f el %.8f diff %.8f\n", d[0],
					jupiter_obs.azimuth*180.0/M_PI, az, az - jupiter_obs.azimuth*180.0/M_PI,
					jupiter_obs.elevation*180.0/M_PI, el, el - jupiter_obs.elevation*180.0/M_PI);

			retval = -1;
		}

		// Increment data line number
		++line;
	}

	return retval;
}
