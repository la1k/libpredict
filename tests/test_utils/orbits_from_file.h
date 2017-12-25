#include <vector>
#include <predict/predict.h>
#include <string>

/**
 * Satellite structure with name and orbital elements.
 **/
typedef struct {
	///Satellite name
	std::string name;
	///Orbital elements
	predict_orbital_elements_t *elements;
	///TLE line 1
	std::string tle_line_1;
	///TLE line 2
	std::string tle_line_2;
} satellite_t;

/**
 * Read TLE file.
 *
 * \param tle_file Filename
 * \return List of orbital elements contained in file
 **/
std::vector<satellite_t> orbital_elements_from_file(const char *tle_file);
