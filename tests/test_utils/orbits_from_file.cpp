#include "orbits_from_file.h"
#include <sstream>

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

