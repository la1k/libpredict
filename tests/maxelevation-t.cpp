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
	//check arguments
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " <testfiles>" << endl;
		return -1;
	}

	//test all provided test files
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

#include <float.h>

//precision of numerical tests
const double EPSILON = 2*FLT_EPSILON;

/**
 * Let predict_at_max_elevation go through a couple of sanity tests:
 * - Check that max elevation time is within AOS/LOS of the pass
 * - Check if the derivative crosses the zero boundary
 * - Check that the elevation is non-zero
 * - Check that rough sampling of the elevation throughout the pass does not yield higher elevations
 * - Check that when the max elevation time is used as the input, the same max elevation time is returned
 * - Check that the same max elevation time is returned regardless of the input time throughout the pass
 * This is done for NUM_PASSES different passes, found from the input time.
 *
 * \param start_time Input start time
 * \param observer Observer
 * \param orbital_elements Orbital elements of satellite
 * \return 0 on success, -1 on failure
 **/
int test_max_elevation(double start_time, predict_observer_t *observer, predict_orbital_elements_t *orbital_elements);

int runtest(const char *filename)
{
	//load testcase
	TestCaseReader testcase;
	testcase.loadFromFile(filename);
	if (!(testcase.containsValidData() && (testcase.containsValidQth()) && (testcase.containsValidTLE()))) {
		fprintf(stderr, "Failed to load testfile: %s\n", filename);
		return -1;
	}

	//get TLE
	char *tle[2];
	testcase.getTLE(tle);

	//create orbit object
	predict_orbital_elements_t *orbital_elements = predict_parse_tle(tle[0], tle[1]);

	if (predict_is_geosynchronous(orbital_elements)) {
		return 0;
	}

	//test at the standard defined QTH
	double start_time = predict_to_julian(testcase.data()[0][0]);
	predict_observer_t *observer = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());
	if (test_max_elevation(start_time, observer, orbital_elements) != 0) {
		fprintf(stderr, "Failed on fixed QTH.\n");
		return -1;
	}

	//test at a QTH that is located exactly on the satellite track, so that elevation should be 90 degrees (and thus potentially a bit problematic with the derivative)
	struct predict_position orbit;
	struct predict_observation obs;
	predict_orbit(orbital_elements, &orbit, start_time);
	observer = predict_create_observer("problematic", orbit.latitude, orbit.longitude, 0);
	double max_elevation_time = start_time;

	//test that the elevation actually is 90.0
	predict_observe_orbit(observer, &orbit, &obs);
	if (!fuzzyCompare(obs.elevation, M_PI/2.0, EPSILON)) {
		fprintf(stderr, "Our QTH does not observe a 90 degrees elevation as it should, something wrong with test: %f\n", obs.elevation*180.0/M_PI);
		return -1;
	}

	//find a time slightly before the pass with elevation less than 0
	double revolution_fraction = 1.0/orbital_elements->mean_motion/12.0;
	while (obs.elevation > 0) {
		start_time -= revolution_fraction;
		predict_orbit(orbital_elements, &orbit, start_time);
		predict_observe_orbit(observer, &orbit, &obs);
	}

	//find AOS/LOS of pass with elevation 90 degrees
	struct predict_observation aos = predict_next_aos(observer, orbital_elements, start_time);
	double next_aos = aos.time;

	struct predict_observation los = predict_next_los(observer, orbital_elements, next_aos);
	double next_los = los.time;
	if (!((max_elevation_time > next_aos) && (max_elevation_time < next_los))) {
		fprintf(stderr, "Error in preparing pass times.\n");
		return -1;
	}

	//test that time predicted from start_time before the 90 degrees pass corresponds to the theoretical maximum elevation
	struct predict_observation max_ele_time_from_start = predict_at_max_elevation(observer, orbital_elements, start_time);
	if (!fuzzyCompare(max_elevation_time, max_ele_time_from_start.time, EPSILON)) {
		fprintf(stderr, "Failed on predicting max elevation corresponding to theoretical 90 degrees time: new = %f, orig = %f\n", max_ele_time_from_start.time, max_elevation_time);
		return -1;
	}

	//test that we get back the theoretical maximum elevation when predicting from theoretical maximum elevation time
	struct predict_observation max_ele_at_max_elevation = predict_at_max_elevation(observer, orbital_elements, max_elevation_time);
	if (!fuzzyCompare(max_elevation_time, max_ele_at_max_elevation.time, EPSILON)) {
		fprintf(stderr, "Failed on predicting the same max elevation as input max elevation for 90 degrees elevation: new = %f, orig = %f\n", max_ele_at_max_elevation.time, max_elevation_time);
		return -1;
	}


	//check that max elevation for times before a pass are correct, and that the rest of the sanity checks are correct for our 90 degrees elevation test qth
	if (test_max_elevation(start_time, observer, orbital_elements) != 0) {
		fprintf(stderr, "Failed on QTH placed directly under the satellite.\n");
		return -1;
	}
	return 0;
}

//threshold for pass length in order to avoid errors in AOS/LOS finding functions triggering an error in the max elevation tests
#define PASS_LENGTH_THRESHOLD (3.0/(24.0*60))

//number of passes to check
#define NUM_PASSES 20

//skip one day between each time point we will predict a pass and check the max elevation
#define TIME_DIFF 1

//number of times we will sample the elevation throughout a pass to roughly check agains the predicted max elevation
#define NUM_TIME_STEPS 10

int test_max_elevation(double start_time, predict_observer_t *observer, predict_orbital_elements_t *orbital_elements)
{
	struct predict_position orbit;
	struct predict_observation obs;

	for (int i=0; i < NUM_PASSES; i++) {
		start_time += i*TIME_DIFF;

		predict_orbit(orbital_elements, &orbit, start_time);
		predict_observe_orbit(observer, &orbit, &obs);

		//one fifth of the total revolution around earth
		double revolution_fraction = 1.0/orbital_elements->mean_motion/5.0;

		//skip to time where current elevation is less than 0 (i.e. skip ahead some time along the revolution around earth)
		if (obs.elevation > 0) {
			struct predict_observation los = predict_next_los(observer, orbital_elements, start_time);
			start_time = los.time + revolution_fraction;
			while (obs.elevation > 0) {
				start_time += revolution_fraction;
				predict_orbit(orbital_elements, &orbit, start_time);
				predict_observe_orbit(observer, &orbit, &obs);
			}
		}

		//get times from next pass
		struct predict_observation aos = predict_next_aos(observer, orbital_elements, start_time);
		predict_julian_date_t next_aos_time = aos.time;
		struct predict_observation los = predict_next_los(observer, orbital_elements, next_aos_time);
		predict_julian_date_t next_los_time = los.time;
		struct predict_observation max_elevation_obs = predict_at_max_elevation(observer, orbital_elements, start_time);
		double max_elevation = max_elevation_obs.elevation;

		//check if max elevation time is within aos/los-times
		if (!((max_elevation_obs.time > next_aos_time) && (max_elevation_obs.time < next_los_time))) {
			fprintf(stderr, "Predicted elevation time not within times for AOS and LOS: predicted time %f, aos %f, los %f\n", max_elevation_obs.time, next_aos_time, next_los_time);
			return -1;
		}

		//check if derivative on either side of the estimated solution crosses the zero boundary
		predict_orbit(orbital_elements, &orbit, max_elevation_obs.time-EPSILON);
		predict_observe_orbit(observer, &orbit, &obs);
		double lower_derivative = obs.elevation_rate;
		predict_orbit(orbital_elements, &orbit, max_elevation_obs.time+EPSILON);
		predict_observe_orbit(observer, &orbit, &obs);
		double upper_derivative = obs.elevation_rate;
		if (!((upper_derivative < 0) && (lower_derivative > 0))) {
			fprintf(stderr, "First derivative of elevation not zero: %f %f %f\n", max_elevation_obs.time, lower_derivative, upper_derivative);
			return -1;
		}

		//check if the found max elevation rate is larger than all other elevations sampled throughout the pass
		for (int i=0; i < NUM_TIME_STEPS; i++) {
			double time = next_aos_time + i*(next_los_time - next_aos_time)/NUM_TIME_STEPS;
			predict_orbit(orbital_elements, &orbit, time);
			predict_observe_orbit(observer, &orbit, &obs);
			if (max_elevation < obs.elevation - EPSILON) {
				fprintf(stderr, "Found elevation through the pass larger than the max elevation: %f, %f\n", max_elevation, obs.elevation);
				return -1;
			}
		}

		//check that the difference between AOS/LOS time is larger than a threshold, and stop the rest of the tests if not (aos/los-functions will detect almost-passes very close to the horizon as actual passes, due to the set thresholds and behavior of the pass finding functions: See issue #18 in the git repository. This will be taken care of among the AOS/LOS tests once the issue is fixed). Predicted max elevation
		//will not fail basic tests above (derivative, actual maximum), but will fail the more nasty tests below since they assume more perfect AOS/LOS timepoints. TODO:
		//Remove this check once issue #18 is fixed, and we have new and more numerically stable aos/los functions.
		if ((next_los_time - next_aos_time) < PASS_LENGTH_THRESHOLD) {
			continue;
		}

		//check if elevation at max elevation is non-zero
		if (max_elevation <= 0) {
			fprintf(stderr, "Predicted max elevation is negative: %f\n", max_elevation);
			return -1;
		}

		//check if max elevation predicted from the max elevation time is valid and the same
		struct predict_observation max_ele_at_max_ele = predict_at_max_elevation(observer, orbital_elements, max_elevation_obs.time);
		if (!((max_ele_at_max_ele.time > next_aos_time) && (max_ele_at_max_ele.time < next_los_time))) {
			fprintf(stderr, "Max elevation predicted at max elevation time not within original pass: %f, %f, %f\n", max_ele_at_max_ele.time, next_aos_time, next_los_time);
			return -1;
		}
		if (!fuzzyCompare(max_elevation_obs.time, max_ele_at_max_ele.time, EPSILON)) {
			fprintf(stderr, "Max elevation time predicted from earlier max elevation time not the same: %f, %f\n", max_elevation_obs.time, max_ele_at_max_ele.time);
			return -1;
		}

		//check if max elevation time predicted from selected time points throughout the pass is consistent
		for (int i=0; i < NUM_TIME_STEPS-1; i++) {
			double time = next_aos_time + i*(next_los_time - next_aos_time)/NUM_TIME_STEPS;
			struct predict_observation new_max_elevation = predict_at_max_elevation(observer, orbital_elements, time);
			if (!fuzzyCompare(new_max_elevation.time, max_elevation_obs.time, EPSILON)) {
				fprintf(stderr, "Failed to predict consistent elevation times: Original time %f, new time %f, at timestep %d\n", max_elevation_obs.time, new_max_elevation.time, i);
				return -1;
			}
		}
	}

	return 0;
}
