#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <predict/predict.h>
#include <iostream>

typedef struct {
	std::string name;
	predict_orbital_elements_t *orbital_elements;
} satellite_t;

/**
 * Parse a TLE file and return as a list of orbital elements.  Based on
 * tle_db_from_file() in flyby, in turn based on ReadDataFiles() from Predict.
 *
 * \param tle_file Path to TLE file
 * \return Parsed orbital elements and associated satellite names as acquired from TLE file.
 **/
std::vector<satellite_t> orbital_elements_from_file(const char *tle_file);

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

	//set fixed start time in UTC close to the epoch of the TLEs
	setenv("TZ", "GMT", 1);
	tzset();

	struct tm timeval = {0};
	timeval.tm_year = 2015-1900;
	timeval.tm_mon = 8;
	timeval.tm_mday = 26;
	timeval.tm_hour = 18;
	timeval.tm_min = 0;

	predict_julian_date_t start_time = predict_to_julian(mktime(&timeval));
	double time_step = 1.0/(24.0*60.0);

	for (int i=0; i < satellite_list.size(); i++) {
		satellite_t satellite = satellite_list[i];
		predict_julian_date_t end_time = start_time + 1.0/satellite.orbital_elements->mean_motion; //mean_motion is the revolutions per day, so the period of the satellite should be 1/mean_motion plus/minus a bit.
		predict_julian_date_t curr_time = start_time;

		std::ofstream file(("sat_track_" + satellite.name + ".dat").c_str());

		while (curr_time < end_time) {
			struct predict_position orbit;
			predict_orbit(satellite.orbital_elements, &orbit, curr_time);
			file << orbit.position[0] << " " << orbit.position[1] << " " << orbit.position[2] << std::endl;
			curr_time += time_step;
		}
		file.close();

		predict_destroy_orbital_elements(satellite.orbital_elements);
	}
}

#define NUM_CHARS_IN_TLE 80

std::vector<satellite_t> orbital_elements_from_file(const char *tle_file)
{
	std::vector<satellite_t> ret_list;

	FILE *fd = fopen(tle_file,"r");
	if (fd == NULL) {
		return ret_list;
	}

	while (feof(fd) == 0) {
		satellite_t satellite;

		char name[NUM_CHARS_IN_TLE] = {0};
		char line1[NUM_CHARS_IN_TLE] = {0};
		char line2[NUM_CHARS_IN_TLE] = {0};

		//read element set
		if (fgets(name, NUM_CHARS_IN_TLE, fd) == NULL) break;
		if (fgets(line1, NUM_CHARS_IN_TLE, fd) == NULL) break;
		if (fgets(line2, NUM_CHARS_IN_TLE, fd) == NULL) break;

		//parse element set
		predict_orbital_elements_t *temp_elements = predict_parse_tle(line1, line2);
		satellite.orbital_elements = temp_elements;
		satellite.name = std::string(name);

		//trim name
		std::stringstream stream(satellite.name);
		stream >> satellite.name;

		ret_list.push_back(satellite);
	}

	fclose(fd);
	return ret_list;
}
