#include "TestCase.h"

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

using namespace std;

TestCase::TestCase()
{
	m_qth_latitude = NAN;
	m_qth_longitude = NAN;
	m_qth_altitude = 0;
}

TestCase::~TestCase()
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

int TestCase::loadFromFile(const char *filename)
{
	
	ifstream file(filename);

	if (!file) {
		cout << "Unable to open file " << filename << endl;
		return -1;
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
				}
			}break;

			case DATA: {
				// Tokenize on =
				vector<string> tokens = tokenize(line, ",");

				// Check size
				if (tokens.size() < 6) continue;
	
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
	
	if (m_data.size() == 0) return -1;
	if (m_tle[0].size() == 0) return -1;
	if (m_tle[1].size() == 0) return -1;
	if (::isnan(m_qth_latitude)) return -1;
	if (::isnan(m_qth_longitude)) return -1;
	
	return 0;
}
	
void TestCase::getTLE(char *tle[2])
{
	tle[0] = new char[m_tle[0].size()];
	tle[1] = new char[m_tle[1].size()];
	strcpy(tle[0], m_tle[0].c_str());
	strcpy(tle[1], m_tle[1].c_str());
}
