#ifndef TEST_CASE_H_
#define TEST_CASE_H_

#include <string>
#include <vector>
#include <limits>

class TestCaseReader
{
public:
	TestCaseReader();
	~TestCaseReader();

	void loadFromFile(const char *filename);

	void getTLE(char *tle[2]);

	double latitude() const {return m_qth_latitude;}
	double longitude() const {return m_qth_longitude;}
	double altitude() const {return m_qth_altitude;}

	double alon() const {return m_alon;};
	double alat() const {return m_alat;};

	std::vector<std::vector<double> > &data() {return m_data;}

	bool containsValidData();
	bool containsValidQth();
	bool containsValidTLE();
	bool containsValidAlonAlat();

private:
	std::string m_tle[2];

	double m_alat;
	double m_alon;

	double m_qth_latitude;
	double m_qth_longitude;
	double m_qth_altitude;

	std::vector<std::vector<double> > m_data;
};

bool fuzzyCompare(const double &x, const double &y, const double &epsilon = std::numeric_limits<double>::epsilon());

/**
 * Check whether input value lies within supplied boundary values. Order of boundary values is irrelevant.
 * \param boundary_value_1 Boundary value 1
 * \param boundary_value_2 Boundary value 2
 * \param compared_value Compared value
 * \param offset Offset added to each boundary value in order to permit larger deviations from the allowed range. Defaults to 0.05, since predict usually is accurate to two decimals
 * \return True if the compared value lies within the boundary values
 **/
bool fuzzyCompareWithBoundaries(const double &boundary_value_1, const double &boundary_value_2, const double &compared_value, double offset = 0.05);

#endif
