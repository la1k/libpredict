#include <stdio.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <predict/predict.h>
#include "testcase_reader.h"
#include <iostream>
extern "C" {
#include "pass_utils.h"
}
using namespace std;

int runtest(const char *filename);

int main(int argc, char **argv)
{
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " <testfiles>" << endl;
		return -1;
	}

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

const double DAYS_PER_SECOND = 1.0/(24.0*60.0*60.0);

bool check_is_root(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, predict_julian_date_t root_candidate)
{
	const double PRECISION = DAYS_PER_SECOND;
	predict_julian_date_t lower = root_candidate - PRECISION;
	predict_julian_date_t upper = root_candidate + PRECISION;

	struct predict_position position;
	struct predict_observation lower_observation, upper_observation;
	predict_orbit(orbital_elements, &position, lower);
	predict_observe_orbit(observer, &position, &lower_observation);

	predict_orbit(orbital_elements, &position, upper);
	predict_observe_orbit(observer, &position, &upper_observation);

	if (lower_observation.elevation * upper_observation.elevation > 0) {
		fprintf(stderr, "Root %f is not valid within %f seconds precision: yielded elevations %f and %f\n", root_candidate, 24*60*60*PRECISION, lower_observation.elevation*180.0/M_PI, upper_observation.elevation*180/M_PI);
		return false;
	}
	return true;
}

int elevation_rate_sign(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, predict_julian_date_t root_candidate)
{
	struct predict_position position;
	struct predict_observation observation;
	predict_orbit(orbital_elements, &position, root_candidate);
	predict_observe_orbit(observer, &position, &observation);
	return observation.elevation_rate/fabs(observation.elevation_rate);
}

bool check_is_aos(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, predict_julian_date_t root_candidate)
{
	if ((elevation_rate_sign(orbital_elements, observer, root_candidate) > 0) && check_is_root(orbital_elements, observer, root_candidate)) {
		return true;
	} else {
		fprintf(stderr, "Input time is not a valid AOS.\n");
		return false;
	}
}

bool check_is_los(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, predict_julian_date_t root_candidate)
{
	if ((elevation_rate_sign(orbital_elements, observer, root_candidate) < 0) && check_is_root(orbital_elements, observer, root_candidate)) {
		return true;
	} else {
		fprintf(stderr, "Input time is not a valid LOS.\n");
		return false;
	}
}

int runtest(const char *filename)
{
	//load testcase data
	TestCaseReader testcase;
	testcase.loadFromFile(filename);
	if (!(testcase.containsValidData() && (testcase.containsValidQth()) && (testcase.containsValidTLE()))) {
		fprintf(stderr, "Test file did not contain valid data: %s\n", filename);
		return -1;
	}

	//prepare TLE and observer
	char *tle[2];
	testcase.getTLE(tle);
	predict_orbital_elements_t *orbital_elements = predict_parse_tle(tle[0], tle[1]);
	predict_observer_t *observer = predict_create_observer("test", testcase.latitude()*M_PI/180.0, testcase.longitude()*M_PI/180.0, testcase.altitude());

	predict_julian_date_t prev_los_time = 0;

	//use the AOS/LOS times in the testcase data files to test the pass stepping methods
	for (unsigned int line_number = 0; line_number < testcase.data().size(); line_number++) {
		std::vector<double> line = testcase.data()[line_number];

		predict_julian_date_t aos_time = line[0];
		predict_julian_date_t los_time = line[1];
		predict_julian_date_t before_pass;

		if (line_number == 0) {
			before_pass = aos_time - DAYS_PER_SECOND;
		} else {
			before_pass = prev_los_time + DAYS_PER_SECOND;
		}

		//check that time before current pass is a valid non-pass timepoint
		struct predict_position orbit;
		struct predict_observation before_pass_observation;
		predict_orbit(orbital_elements, &orbit, before_pass);
		predict_observe_orbit(observer, &orbit, &before_pass_observation);

		if (before_pass_observation.elevation >= 0) {
			fprintf(stderr, "Observation before pass is not valid, elevation %f. Data line %d.\n", before_pass_observation.elevation, line_number);
			return -1;
		}

		//check that testdata AOS and LOS times are valid according to the requirements
		if (!check_is_aos(orbital_elements, observer, aos_time) || !check_is_los(orbital_elements, observer, los_time)) {
			fprintf(stderr, "Failed testcase data check at line %d: testcase data is not valid\n", line_number);
			return -1;
		}

		/*
		const int NUM_STARTING_POINTS = 10;
		for (int i=0; i < NUM_STARTING_POINTS; i++) {
			predict_julian_date_t start_time = i*(aos_time - before_pass)/NUM_STARTING_POINTS + before_pass;

			//try to step into the pass
			predict_julian_date_t pass_time = step_pass(observer, orbital_elements, start_time, POSITIVE_DIRECTION, STEP_UNTIL_PASS);
			if (pass_time < aos_time) {
				fprintf(stderr, "Did not step into the pass at line %d, starting point %d\n", line_number, i);
				return -1;
			}
			if (pass_time > los_time) {
				fprintf(stderr, "Stepped over the pass at line %d, starting point %d\n", line_number, i);
				return -1;
			}

			//try to step out of the pass
			if (line_number < testcase.data().size()) {
				std::vector<double> next_line = testcase.data()[line_number+1];
				predict_julian_date_t start_time = i*(los_time - aos_time)/NUM_STARTING_POINTS + aos_time;

				predict_julian_date_t out_of_pass_time = step_pass(observer, orbital_elements, start_time, POSITIVE_DIRECTION, STEP_OUT_OF_PASS);

				predict_julian_date_t next_aos_time = predict_to_julian(next_line[0]);
				if (out_of_pass_time < los_time) {
					fprintf(stderr, "Did not step out of the pass at line %d, starting point %d\n", line_number, i);
					return -1;
				}

				if (out_of_pass_time > next_aos_time) {
					fprintf(stderr, "Stepped out of the pass, but stepped over the next pass. Line %d, starting point %d\n", line_number, i);
					return -1;
				}
			}
		}*/

		line_number++;
		prev_los_time = los_time;
	}

	return 0;
}
