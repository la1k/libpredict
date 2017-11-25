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
		double time = d[0];
		double az = d[1];
		double el = d[2];
		double dec = d[3];
		double gha = d[4];
		double ra = d[5];

		// Compare values within (time - 1, time + 1) (i.e. up time + 1, but not including time + 1)
		// (since we don't know the exact time predict generated its data, only within an error of 1 second)
		const int DIFF = 1;
		struct predict_observation moon_obs_lower;
		struct predict_observation moon_obs_upper;

		// Lower bound
		predict_observe_moon(obs, predict_to_julian(time), &moon_obs_lower);

		// Upper bound
		predict_observe_moon(obs, predict_to_julian(time + DIFF), &moon_obs_upper);

		// Check values
		string failed = "";
		if (!fuzzyCompareWithBoundaries(moon_obs_lower.azimuth*180.0/M_PI, moon_obs_upper.azimuth*180.0/M_PI, az)) {
			failed += "(azimuth)";
		}
		if (!fuzzyCompareWithBoundaries(moon_obs_lower.elevation*180.0/M_PI, moon_obs_upper.elevation*180.0/M_PI, el)) {
			failed += "(elevation)";
		}

		//calculate RA, dec and GHA
		double dec_lower = predict_moon_declination(predict_to_julian(time))*180.0/M_PI;
		double dec_upper = predict_moon_declination(predict_to_julian(time + DIFF))*180.0/M_PI;
		if (!fuzzyCompareWithBoundaries(dec_lower, dec_upper, dec)) {
			failed += "(declination)";
		}

		double ra_lower = predict_moon_ra(predict_to_julian(time))*180.0/M_PI;
		double ra_upper = predict_moon_ra(predict_to_julian(time + DIFF))*180.0/M_PI;
		if (!fuzzyCompareWithBoundaries(ra_lower, ra_upper, ra)) {
			failed += "(right ascension)";
		}

		double gha_lower = predict_moon_gha(predict_to_julian(time))*180.0/M_PI;
		double gha_upper = predict_moon_gha(predict_to_julian(time + DIFF))*180.0/M_PI;
		if (!fuzzyCompareWithBoundaries(gha_lower, gha_upper, gha)) {
			failed += "(GHA)";
		}

		// Failed?
		if (failed != "") {
			cout << filename << ": failed at data line " << line << ": " << failed << endl;

			printf("%.8f, %.8f/%.8f/%.8f, %.8f/%.8f/%.8f, "
					"%.8f/%.8f/%.8f, "
					"%.8f/%.8f/%.8f, "
					"%.8f/%.8f/%.8f,"
					"\n", time,
					moon_obs_lower.azimuth*180.0/M_PI, az, moon_obs_upper.azimuth*180.0/M_PI,
					moon_obs_lower.elevation*180.0/M_PI, el, moon_obs_upper.elevation*180.0/M_PI,
					dec_lower, dec, dec_upper,
					ra_lower, ra, ra_upper,
					gha_lower, gha, gha_upper
			      );


			retval = -1;
		}

		// Increment data line number
		++line;
	}

	return retval;
}
