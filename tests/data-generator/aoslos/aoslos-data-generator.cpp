#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <predict/predict.h>

#include <test_utils/orbits_from_file.h>

#include <defs.h>
extern "C" {
#include "pass_utils.h"
}

void print_progress_bar(std::string action, long last_progress, long progress, long total) {
	if (progress == 0) {
		fprintf(stderr, "%s\n\n", action.c_str());
	}

	const int MAX_NUM_SIGNS = 20;
	int num_signs = progress*1.0/(total*1.0)*MAX_NUM_SIGNS;
	int prev_num_signs = last_progress*1.0/(total*1.0)*MAX_NUM_SIGNS;

	if ((prev_num_signs != num_signs) || (progress == 0)) {
		char str[MAX_NUM_SIGNS+3] = {0};
		str[0] = '[';
		str[MAX_NUM_SIGNS] = ']';
		for (int i=1; i < MAX_NUM_SIGNS; i++) {
			if (i < num_signs) {
				str[i] = '#';
			} else {
				str[i] = ' ';
			}
		}
		fprintf(stderr, "\033[1A%s %ld of %ld\n", str, progress, total);
	}
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		fprintf(stderr, "Usage: %s <tle-file>\n", argv[0]);
		exit(1);
	}

	//parse TLE file
	std::vector<satellite_t> satellite_list = orbital_elements_from_file(argv[1]);
	if (satellite_list.size() <= 0) {
		fprintf(stderr, "No TLEs found in input file.\n");
		exit(1);
	}

	//create observers
	std::vector<predict_observer_t*> observers = {
						      predict_create_observer("QTH_1", 63.4375*M_PI/180.0, 10.37500*M_PI/180.0, 0),
						      predict_create_observer("QTH_2", -63.4375*M_PI/180.0, 10.37500*M_PI/180.0, 0),
						      predict_create_observer("QTH_3", 63.4375*M_PI/180.0, -10.37500*M_PI/180.0, 0),
						      predict_create_observer("QTH_4", -63.4375*M_PI/180.0, -10.37500*M_PI/180.0, 0),
						     };

	//create test data for each satellite and observer combination
	for (int i=0; i < satellite_list.size(); i++) {
		satellite_t satellite = satellite_list[i];
		if (predict_is_geosynchronous(satellite.elements)) {
			continue;
		}

		double tle_epoch_year = satellite.elements->epoch_year + 2000.0;
		double tle_epoch_day = satellite.elements->epoch_day;

		//start predictions at one day after the epoch day
		struct tm timeval = {0};
		timeval.tm_year = tle_epoch_year-1900;
		timeval.tm_yday = tle_epoch_day+1;
		time_t start_unix_time = mktime(&timeval);
		predict_julian_date_t start_time = predict_to_julian(start_unix_time);

		//should predict to to ten times whichever is the largest of earth's rotational period and the satellite's mean orbital period
		double orbital_period = 1.0/satellite.elements->mean_motion;
		double earth_period = 1.0/EARTH_ROTATIONS_PER_SIDERIAL_DAY;
		double canonical_period = (orbital_period < earth_period ? orbital_period : earth_period);
		const double NUM_PERIODS = 10;
		predict_julian_date_t end_time = start_time + NUM_PERIODS*canonical_period;
		for (int j=0; j < observers.size(); j++) {
			predict_observer_t *observer = observers[j];
			if (!predict_aos_happens(satellite.elements, observer->latitude)) {
				continue;
			}

			//loop through elevation curves at a fine time grid, collect times for AOS and LOS
			double time_step = canonical_period/(24.0*60.0);
			predict_julian_date_t time = start_time;
			struct predict_position orbit;
			struct predict_observation observation;

			predict_orbit(satellite.elements, &orbit, time);
			predict_observe_orbit(observer, &orbit, &observation);
			struct predict_observation prev = observation;

			bool in_beginning_pass = (prev.elevation >= 0);
			double non_pass_time_before_pass = start_time;

			std::vector<double> aos_timepoints;
			std::vector<double> los_timepoints;

			long tot_iterations = (end_time - start_time)/time_step;
			long iterations = 0;

			std::ofstream elevation_file("elevations_" + satellite.name + "_qth_" + std::string(observer->name) + ".dat");
			std::ofstream event_file("events_" + satellite.name + "_qth_" + std::string(observer->name) + ".dat");

			while ((time <= end_time) || (aos_timepoints.size() != los_timepoints.size())) {
				predict_orbit(satellite.elements, &orbit, time);
				predict_observe_orbit(observer, &orbit, &observation);
				double print_time = time - start_time;

				predict_julian_date_t root_time = -1;

				if (prev.elevation*observation.elevation < 0) {
					//this is a root: fine-tune using bisection method
					root_time = bisection_method(satellite.elements, observer, prev.time, observation.time);
				}

				if ((prev.elevation < 0) && (observation.elevation >= 0) && !in_beginning_pass) {
					aos_timepoints.push_back(root_time);
					event_file << print_time << " " << observation.elevation << std::endl;
				} else if ((prev.elevation >= 0) && (observation.elevation < 0) && !in_beginning_pass) {
					los_timepoints.push_back(root_time);
					event_file << print_time << " " << observation.elevation << std::endl;
				}
				elevation_file << print_time << " " << observation.elevation << std::endl;

				if (in_beginning_pass && observation.elevation < 0) {
					in_beginning_pass = false;
					non_pass_time_before_pass = time;
				}

				prev = observation;

				print_progress_bar(satellite.name + ", " + std::string(observer->name), iterations-1, iterations, tot_iterations);
				time += time_step;
				iterations++;
			}

			elevation_file.close();

			//output as testcase data to file
			std::ofstream output_file("aoslos_" + satellite.name + "_qth_" + std::string(observer->name) + ".test");
			output_file.precision(20);
			output_file << "[tle]" << std::endl
				    << satellite.tle_line_1
				    << satellite.tle_line_2 << std::endl
				    << "[qth]" << std::endl
				    << "lat=" << observer->latitude*180.0/M_PI << std::endl
				    << "lon=" << observer->longitude*180.0/M_PI << std::endl
				    << "alt=" << observer->altitude*180.0/M_PI << std::endl << std::endl
				    << "[data]" << std::endl
				    << non_pass_time_before_pass << std::endl;

			for (int data_ind=0; data_ind < aos_timepoints.size(); data_ind++) {
				output_file << aos_timepoints[data_ind] << "," << los_timepoints[data_ind] << std::endl;
			}

			output_file.close();
		}
	}
}

