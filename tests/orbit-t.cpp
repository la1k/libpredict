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
	predict_orbital_elements_t *orbital_elements = predict_parse_tle(tle[0], tle[1]);

	// Used in lower bound in value check
	struct predict_position orbit_lower;

	// Used in upper bound in value check
	struct predict_position orbit_upper;

	// Create observer object
	predict_observer_t *obs = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		return -1;
	}

	bool check_squint_angle = testcase.containsValidAlonAlat() && (orbital_elements->ephemeris == EPHEMERIS_SDP4);
	
	// Test
	int retval = 0;
	int line = 1;
	for (size_t i=0; i < testcase.data().size(); i++) {
		std::vector<double> d = testcase.data()[i];
		double time = d[0];
		double lat = d[1];
		double lon = d[2];
		double alt = d[3];
		double az = d[4];
		double el = d[5];
		double doppler = d[6];
		double squint = d[7];
		double phase = d[8];
		long revolutions = d[9];
		double footprint = d[10];
		double range = d[11];
		double velocity = d[12];
		double visibility = d[13];
		double eclipse_depth = d[14];

		bool is_in_sunlight = false;
		bool is_visible = false;

		// Parse visibility value
		if (visibility == 0.0) {
			is_in_sunlight = false;
			is_visible = false;
		} else if (visibility == 1.0) {
			is_in_sunlight = true;
			is_visible = false;
		} else if (visibility == 2.0) {
			is_in_sunlight = true;
			is_visible = true;
		}

		// Compare values within (time - 1, time + 1) (i.e. up time + 1, but not including time + 1)
		// (since we don't know the exact time predict generated its data, only within an error of 1 second)
		const int DIFF = 1;
		struct predict_observation orbit_obs_lower;
		struct predict_observation orbit_obs_upper;

		// Lower bound
		predict_orbit(orbital_elements, &orbit_lower, predict_to_julian(time));
		predict_observe_orbit(obs, &orbit_lower, &orbit_obs_lower);

		// Upper bound
		predict_orbit(orbital_elements, &orbit_upper, predict_to_julian(time + DIFF));
		predict_observe_orbit(obs, &orbit_upper, &orbit_obs_upper);

		// Check values
		string failed = "";

		// Lat, lon, alt
		if (!fuzzyCompareWithBoundaries(orbit_lower.latitude*180.0/M_PI, orbit_upper.latitude*180/M_PI, lat)) {
			failed += "(latitude)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_lower.longitude*180.0/M_PI, orbit_upper.longitude*180/M_PI, lon)) {
			failed += "(longitude)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_lower.altitude, orbit_upper.altitude, alt)) {
			failed += "(altitude)";
		}

		// Azi, ele
		if (!fuzzyCompareWithBoundaries(orbit_obs_lower.azimuth*180.0/M_PI, orbit_obs_upper.azimuth*180.0/M_PI, az)) {
			failed += "(azimuth)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_obs_lower.elevation*180.0/M_PI, orbit_obs_upper.elevation*180.0/M_PI, el)) {
			failed += "(elevation)";
		}

		// Doppler shift, footprint, range, velocity
		double frequency = 100.0e06;  //predict outputs a weird factor instead of the actual doppler shift (since sign depends on whether it is uplink or downlink frequency), so can set this fixed frequency in order to get the same factor from libpredict.
        // uplink and downlink correction is different as we want to anticipate the doppler observed by the satellite such that the satellite receives the uplink signal at nominal frequency
		double doppler_lower = predict_doppler_shift(&orbit_obs_lower, frequency);
		double doppler_upper = predict_doppler_shift(&orbit_obs_upper, frequency);
		if (!fuzzyCompareWithBoundaries(doppler_lower, doppler_upper, doppler)) {
			failed += "(doppler)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_lower.footprint, orbit_upper.footprint, footprint, 1)) {
			failed += "(footprint)";
		}
		if (!fuzzyCompareWithBoundaries(orbit_obs_lower.range, orbit_obs_upper.range, range)) {
			failed += "(range)";
		}
		double velocity_lower = sqrt(pow(orbit_lower.velocity[0], 2.0) + pow(orbit_lower.velocity[1], 2.0) + pow(orbit_lower.velocity[2], 2.0));
		double velocity_upper = sqrt(pow(orbit_upper.velocity[0], 2.0) + pow(orbit_upper.velocity[1], 2.0) + pow(orbit_upper.velocity[2], 2.0));
		if (!fuzzyCompareWithBoundaries(velocity_lower, velocity_upper, velocity)) {
			failed += "(velocity)";
		}

		// Eclipse depth
		if (!fuzzyCompareWithBoundaries(orbit_lower.eclipse_depth*180.0/M_PI, orbit_upper.eclipse_depth*180.0/M_PI, eclipse_depth)) {
			failed += "(eclipse_depth)";
		}

		// Visibility status
		if (!(orbit_lower.eclipsed) != is_in_sunlight) {
			failed += "(eclipsed)";
		}
		if (orbit_obs_lower.visible != is_visible) {
			failed += "(visibility)";
		}

		// Squint angle
		double squint_angle_lower, squint_angle_upper;
		if (check_squint_angle) {
			squint_angle_lower = predict_squint_angle(obs, &orbit_lower, testcase.alon()*M_PI/180.0, testcase.alat()*M_PI/180.0)*180.0/M_PI;
			squint_angle_upper = predict_squint_angle(obs, &orbit_upper, testcase.alon()*M_PI/180.0, testcase.alat()*M_PI/180.0)*180/M_PI;
			if (!fuzzyCompareWithBoundaries(squint_angle_lower, squint_angle_upper, squint)) {
				failed += "(squint)";
			}
		}

		// Phase
		if (!fuzzyCompareWithBoundaries(orbit_lower.phase*180.0/M_PI, orbit_upper.phase*180.0/M_PI, phase)) {
			failed += "(phase)";
		}

		// Revolutions
		if (!fuzzyCompareWithBoundaries(orbit_lower.revolutions, orbit_upper.revolutions, revolutions)) {
			failed += "(revolutions)";
		}

		// Failed?
		if (failed != "") {
			cout << filename << ": failed at data line " << line << ": " << failed << endl;

			printf("%.8f, %.8f/%.8f/%.8f, %.8f/%.8f/%.8f, %.3f/%.3f/%.3f, %.3f/%.3f/%.3f, %.3f/%.3f/%.3f", time,
					orbit_lower.latitude*180.0/M_PI, lat, orbit_upper.latitude*180.0/M_PI,
					orbit_lower.longitude*180.0/M_PI, lon, orbit_upper.longitude*180.0/M_PI,
					orbit_lower.altitude, alt, orbit_upper.altitude,
					orbit_obs_lower.azimuth*180.0/M_PI, az, orbit_obs_upper.azimuth*180.0/M_PI,
					orbit_obs_lower.elevation*180.0/M_PI, el, orbit_obs_upper.elevation*180.0/M_PI);
			printf(", %.8f/%.8f/%.8f", doppler_lower, doppler, doppler_upper);
			printf(", %.8f/%.8f/%.8f", orbit_lower.footprint, footprint, orbit_upper.footprint);
			printf(", %.8f/%.8f/%.8f", orbit_obs_lower.range, range, orbit_obs_upper.range);
			printf(", %.8f/%.8f/%.8f", velocity_lower, velocity, velocity_upper);
			printf(", %.8f/%.8f/%.8f", orbit_lower.eclipse_depth*180.0/M_PI, eclipse_depth, orbit_upper.eclipse_depth*180.0/M_PI);
			printf(", %d/%d", !(orbit_lower.eclipsed), is_in_sunlight);
			printf(", %d/%d", orbit_obs_lower.visible, is_visible);

			if (check_squint_angle) {
				printf(", %.8f/%.8f/%.8f", squint_angle_lower, squint, squint_angle_upper);
			}
			printf(", %.8f/%.8f/%.8f, %ld/%ld/%ld", orbit_lower.phase*180.0/M_PI, phase, orbit_upper.phase*180.0/M_PI, orbit_lower.revolutions, revolutions, orbit_upper.revolutions);
			printf("\n");

			retval = -1;
		}

		// Increment data line number
		++line;
	}

	return retval;
}
