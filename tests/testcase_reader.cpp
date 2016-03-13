#include "testcase_reader.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <string.h>

#include <cmath>

using namespace std;

TestCaseReader::TestCaseReader()
{
	m_qth_latitude = NAN;
	m_qth_longitude = NAN;
	m_qth_altitude = 0;

	m_alat = NAN;
	m_alon = NAN;
}

TestCaseReader::~TestCaseReader()
{

}
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

vector<string> tokenize(const string& str,const string& delimiters)
{
	vector<string> tokens;
	string::size_type delimPos = 0, tokenPos = 0, pos = 0;

	if(str.length()<1)  return tokens;
	while(1){
		delimPos = str.find_first_of(delimiters, pos);
		tokenPos = str.find_first_not_of(delimiters, pos);
		if(string::npos != delimPos){
			if(string::npos != tokenPos){
				if(tokenPos<delimPos){
					tokens.push_back(str.substr(pos,delimPos-pos));
				}else{
					tokens.push_back("");
				}
			}else{
				tokens.push_back("");
			}
			pos = delimPos+1;
		}else{
			if(string::npos != tokenPos){
				tokens.push_back(str.substr(pos));
			}else{
				tokens.push_back("");	
			}
			break;
		}
	}
	return tokens;
}

void TestCaseReader::loadFromFile(const char *filename)
{
	
	ifstream file(filename);

	if (!file) {
		cout << "Unable to open file " << filename << endl;
		return;
	}
	
	enum State {
		SEARCH,
		TLE,
		QTH,
		DATA,
	} state = SEARCH;

	int count = 0;

	string line;
	while (getline(file, line)) {

		// Trim and ignore # and empty lines
		line = trim(line);
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;
		
		// Detect block change
		if (line == "[tle]") {
			state = TLE;
			count = 0;
			continue;
		} else if (line == "[qth]") {
			state = QTH;
			count = 0;
			continue;
		} else if (line == "[data]") {
			state = DATA;
			count = 0;
			continue;
		}
		
		// Handle state
		switch (state) {

			case TLE:
				if (count <= 1) {
					m_tle[count++] = line;
				} else {
					state = SEARCH;
					count = 0;
				}
			break;

			case QTH: {
				// Tokenize on =
				vector<string> tokens = tokenize(line, "=");

				// Check size
				if (tokens.size() != 2) continue;

				// Trim key and value
				string key = trim(tokens[0]);
				string value = trim(tokens[1]);

				// Skip if empty
				if (key.size() == 0) continue;
				if (value.size() == 0) continue;
				
				// Handle keywords
				if (key == "lat" || key == "latitude") {
					stringstream(value) >> m_qth_latitude;

				} else if (key == "lon" || key == "longitude") {
					stringstream(value) >> m_qth_longitude;
				} else if (key == "alt" || key == "altitude") {
					stringstream(value) >> m_qth_altitude;
				} else if (key == "alon"){
					if (value != "No alon") {
						stringstream(value) >> m_alon;
					}
				} else if (key == "alat"){
					if (value != "No alat") {
						stringstream(value) >> m_alat;
					}
				}
			} break;

			case DATA: {
				// Tokenize on =
				vector<string> tokens = tokenize(line, ",");

				vector<double> d;
				// Convert to double vector
				for (int i=0;i<(int)tokens.size();i++) {
					// Trim token spaces
					tokens[i] = trim(tokens[i]);
					
					// Convert to double
					double val;
					stringstream(tokens[i]) >> val;

					// Append
					d.push_back(val);
				}

				// Save in data vector
				m_data.push_back(d);
				
			} break;
			default: break;
		}
	}

	file.close();
}
	
void TestCaseReader::getTLE(char *tle[2])
{
	tle[0] = new char[m_tle[0].size()];
	tle[1] = new char[m_tle[1].size()];
	strcpy(tle[0], m_tle[0].c_str());
	strcpy(tle[1], m_tle[1].c_str());
}

bool TestCaseReader::containsValidData()
{
	return (m_data.size() != 0);
}

bool TestCaseReader::containsValidQth()
{
	return (!(std::isnan(m_qth_latitude)) && !(std::isnan(m_qth_longitude)));
}

bool TestCaseReader::containsValidTLE()
{
	return ((m_tle[0].size() != 0) && (m_tle[1].size() != 0));
}

bool TestCaseReader::containsValidAlonAlat()
{
	return (!(std::isnan(m_alon)) && !(std::isnan(m_alat)));
}

bool fuzzyCompare(const double &x, const double &y, const double &epsilon)
{
	return fabs(x - y) < epsilon;
}

bool fuzzyCompareWithBoundaries(const double &input_value_1, const double &input_value_2, const double &compared_value, double offset)
{
	double lower, upper;
	if (input_value_2 > input_value_1)
	{
		lower = input_value_1;
		upper = input_value_2;
	}
	else
	{
		lower = input_value_2;
		upper = input_value_1;
	}
	return (compared_value < upper + offset) && (compared_value > lower - offset);
}
