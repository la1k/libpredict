#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <predict/predict.h>

#include "testcase_reader.h"

#include <sstream>
#include <iostream>
#include <fstream>

/**
 * Satellite structure with name and orbital elements.
 **/
typedef struct {
	///Satellite name
	std::string name;
	///Orbital elements
	predict_orbital_elements_t *elements;
} satellite_t;

/**
 * Read TLE file.
 *
 * \param tle_file Filename
 * \return List of orbital elements contained in file
 **/
std::vector<satellite_t> orbital_elements_from_file(const char *tle_file);

/**
 * Read list of satellite numbers from file.
 *
 * \param filename Filename
 * \return Satellite number list, where retmap[satnum] returns true if satellite number is present in the file
 **/
std::map<long, bool> satcat_list_from_file(const char *filename);

/**
 * Print "geostationary" or "non-geostationary" depending on input. Used for printing convenience in error messages in main().
 *
 * \param geos Geostationary/non-geostationary
 * \return "geostationary" if geos is true, "non-geostationary" otherwise
 **/
const char* geos_string(bool geos);

int main(int argc, char **argv)
{
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " tle-file [optional: geostationary-satcats]" << std::endl;
		return 1;
	}

	std::vector<satellite_t> orbital_elements = orbital_elements_from_file(argv[1]);
	std::map<long, bool> geostationary;

	//If satellite number list is not specified, we assume that all satellites in TLE file should be geostationary
	bool all_geostationary = false;
	if (argc < 3) {
		all_geostationary = true;
	} else {
		geostationary = satcat_list_from_file(argv[2]);
	}
	bool failed = false;

	for (size_t i=0; i < orbital_elements.size(); i++) {
		predict_orbital_elements_t *elements = orbital_elements[i].elements;
		std::string name = orbital_elements[i].name;
		bool predict_geostationary = predict_is_geosynchronous(elements);
		bool is_geostationary = true;
		if (!all_geostationary) {
			//check against list over satellites that should be geostationary
			is_geostationary = geostationary[elements->satellite_number];
		}

		if (predict_geostationary != is_geostationary) {
			fprintf(stderr, "Predicted %s, but should be %s ", geos_string(predict_geostationary), geos_string(is_geostationary));
			fprintf(stderr, "for %s (%i): meanmo=%f, ecc=%f, inclin=%f\n", name.c_str(), elements->satellite_number, elements->mean_motion, elements->eccentricity, elements->inclination);
			failed = true;

		}
	}
	return failed;
}

#define NUM_CHARS_IN_TLE 80
std::vector<satellite_t> orbital_elements_from_file(const char *tle_file)
{
	std::vector<satellite_t> ret_list;

	FILE *fd = fopen(tle_file, "r");
	if (fd == NULL) {
		fprintf(stderr, "Failed to open TLE file.\n");
		exit(1);
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
		satellite.elements = temp_elements;
		satellite.name = std::string(name);

		//trim name
		std::stringstream stream(satellite.name);
		stream >> satellite.name;

		ret_list.push_back(satellite);
	}

	if (ret_list.size() == 0) {
		fprintf(stderr, "Failed to read TLEs from TLE file.\n");
		exit(1);
	}

	fclose(fd);
	return ret_list;
}

std::map<long, bool> satcat_list_from_file(const char *filename)
{
	std::ifstream file;
	file.open(filename);

	if (file.fail()) {
		fprintf(stderr, "Failed to open satellite numbers file.\n");
		exit(1);
	}

	std::map<long, bool> ret_list;
	while (!file.eof()) {
		long sat_cat;
		file >> sat_cat;
		ret_list[sat_cat] = true;
	}
	return ret_list;
}

const char* geos_string(bool geos)
{
	if (geos) {
		return "geostationary";
	} else {
		return "non-geostationary";
	}
}
