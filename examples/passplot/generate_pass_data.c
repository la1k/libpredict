#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include <predict/predict.h>

double rad_to_deg(double rad)
{
	double deg = rad*180.0/M_PI;
	return deg;
}

int main(int argc, char **argv)
{
	//TLE of OSCAR 7 obtained at 2017-02-27
	const char *oscar7_tle_line_1 = "1 07530U 74089B   17058.02491442 -.00000048  00000-0 -22049-4 0  9995";
	const char *oscar7_tle_line_2 = "2 07530 101.6163  28.9438 0011984 174.4353 227.0960 12.53625643935054";
	predict_orbital_elements_t *orbital_elements = predict_parse_tle(oscar7_tle_line_1, oscar7_tle_line_2);
	float latitude = 63.9;
	float longitude = 10.9;
	predict_observer_t *observer = predict_create_observer("LA1K", latitude*M_PI/180.0, longitude*M_PI/180.0, 10);

	//force UTC time in order to fix correct time below (in general not necessary for correct libpredict operation, as libpredict does conversions to and from timezones automatically)
	setenv("TZ", "GMT", 1);
	tzset();

	//fixed time close to TLE epoch date so that the TLE is valid.
	//the fixed time corresponds to 2017-02-28 11:30 UTC
	struct tm timeval = {0};
	timeval.tm_year = 2017-1900;
	timeval.tm_mon = 1; //this is february
	timeval.tm_mday = 28;
	timeval.tm_hour = 11;
	timeval.tm_min = 30;

	//predict next AOS and LOS after defined time
	predict_julian_date_t start_time = predict_to_julian(mktime(&timeval));
	struct predict_observation aos, los;
	aos = predict_next_aos(observer, orbital_elements, start_time);
	los = predict_next_los(observer, orbital_elements, start_time);

	//calculate properties between AOS and LOS
	struct predict_position orbit;
	struct predict_observation observation;

	double timestep = 1.0/(24.0*60.0*60);
	predict_julian_date_t curr_time = aos.time - 30*60*timestep;
	predict_julian_date_t end_time = los.time + 30*60*timestep;

	double freq = 1000;

	while (curr_time < end_time) {
		predict_orbit(orbital_elements, &orbit, curr_time);
		predict_observe_orbit(observer, &orbit, &observation);
		double shift = predict_doppler_shift(&observation, freq);

		//print independent properties
		printf("%f %f ", rad_to_deg(orbit.latitude), rad_to_deg(orbit.longitude));

		//print observed properties
		printf("%f %f %f ", rad_to_deg(observation.elevation), rad_to_deg(observation.azimuth), freq + shift);
		printf("%f\n", rad_to_deg(observation.elevation_rate));
		curr_time += timestep;
	}

	predict_destroy_observer(observer);
	predict_destroy_orbital_elements(orbital_elements);
}
