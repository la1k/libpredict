#ifndef TEST_CASE_H_
#define TEST_CASE_H_

#include <string>
#include <vector>

class TestCase
{
public:
	TestCase();
	~TestCase();

	int loadFromFile(const char *filename);

	void getTLE(char *tle[2]);

	double latitude() const {return m_qth_latitude;}
	double longitude() const {return m_qth_longitude;}
	double altitude() const {return m_qth_altitude;}

	std::vector<std::vector<double> > &data() {return m_data;}

private:
	std::string m_tle[2];

	double m_qth_latitude;
	double m_qth_longitude;
	double m_qth_altitude;

	std::vector<std::vector<double> > m_data;

};

#endif
