#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <predict/predict.h>

double observer_next_sunset(const predict_observer_t *observer, double time, struct predict_observation *obs)
{
	struct predict_observation sun;

	// Scan for first elevation > 0 (t0)
	double t0 = time;
	predict_observe_sun(observer, t0, &sun);
	while (sun.elevation < 0) {
		t0 += 1.0/24.0/3600.0;
		predict_observe_sun(observer, t0, &sun);
	}

	//Scan for elevation < 0 (t1);
	double t1 = t0 + 1.0/24.0/3600.0;
	predict_observe_sun(observer, t1, &sun);
	while (sun.elevation > 0) {
		t1 += 1.0/24.0/3600.0;
		predict_observe_sun(observer, t1, &sun);
	}

	while ( fabs(t0-t1) > 0.01/24/3600 ) {

		time = (t0 + t1) / 2.0;
		predict_observe_sun(observer, time, &sun);

		if (sun.elevation > 0) t0 = time;
		else t1 = time;

	}

	if (obs != NULL) {
		memcpy(obs, &sun, sizeof(struct predict_observation));
	}

	return time;
}

double observer_next_sunrise(const predict_observer_t *observer, double time, struct predict_observation *obs)
{
	struct predict_observation sun;

	// Scan for first elevation < 0 (t0)
	double t0 = time;
	predict_observe_sun(observer, t0, &sun);
	while (sun.elevation > 0) {
		t0 += 1.0/24.0/3600.0;
		predict_observe_sun(observer, t0, &sun);
	}

	//Scan for elevation > 0 (t1);
	double t1 = t0 + 1.0/24.0/3600.0;
	predict_observe_sun(observer, t1, &sun);
	while (sun.elevation < 0) {
		t1 += 1.0/24.0/3600.0;
		predict_observe_sun(observer, t1, &sun);
	}

	while ( fabs(t0-t1) > 0.01/24/3600 ) {

		time = (t0 + t1) / 2.0;
		predict_observe_sun(observer, time, &sun);

		if (sun.elevation < 0) t0 = time;
		else t1 = time;

	}

	if (obs != NULL) {
		memcpy(obs, &sun, sizeof(struct predict_observation));
	}

	return time;
}

int main(int argc, char **argv)
{

	// Create observer object
	predict_observer_t *obs = predict_create_observer("Me", 59.95*M_PI/180.0, 10.75*M_PI/180.0, 0);
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		exit(1);
	}

	predict_julian_date_t curr_time = predict_to_julian(time(NULL));

	struct predict_observation sun;
	double sunset = observer_next_sunset(obs, curr_time, &sun);

	// Convert to hour, minute, seconds
	double timeto = (sunset - curr_time)*24*3600;
	int h = timeto / 3600;
	int m = (timeto-h*3600) / 60;
	int s = ((int)timeto)%60;

	time_t t = (time_t)(86400.0 * (sunset + 3651.0));
	printf("Next sunset in %02i:%02i:%02i, azimuth=%.1f, at UTC %s", h, m, s, sun.azimuth*180.0/M_PI, asctime(gmtime(&t)));

	double sunrise = observer_next_sunrise(obs, curr_time, &sun);

	// Convert to hour, minute, seconds
	timeto = (sunrise - curr_time)*24*3600;
	h = timeto / 3600;
	m = (timeto-h*3600) / 60;
	s = ((int)timeto)%60;

	t = (time_t)(86400.0 * (sunrise + 3651.0));
	printf("Next sunrise in %02i:%02i:%02i, azimuth=%.1f, at UTC %s", h, m, s, sun.azimuth*180.0/M_PI, asctime(gmtime(&t)));

	// Free memory
	predict_destroy_observer(obs);

	return 0;
}
