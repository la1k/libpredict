#!/bin/bash

## Master file for generating test data from predict for use in validating libpredict's API functions. 
## Uses the UDP server functionality in predict for generating data.
## Will compile an UDP client and use this for communicating. 
##
## See list at the bottom of this file for files that are generated. 
##
## Note that the granularity of the time is accurate to a second, while the obtained data can be 
## anywhere from the obtained timestamp up till the next second. predict also uses only two decimals in
## the output. 
##
## Also note that due to internal restrictions in predict, TLE files should be restricted to 24 satellites, as 
## there otherwise could be undefined behavior in the UDP calls to the predict server. 


# Generate QTH lines for testcase file from predict QTH file.
function get_qth_string(){
	qth_file="$1"
	echo "[qth]"
	echo "lat=$(sed '2!d' $qth_file | sed -rn 's/\s//p')"
	lon=$(sed '3!d' $qth_file | sed -rn 's/\s//p')
	
	#convert longitude from N/W to N/E
	lon=$(echo "-1*$lon" | bc)

	echo "lon=$lon"
	echo "alt=$(sed '4!d' $qth_file | sed -rn 's/\s//p')"
}

# Generate testcase data for a given satellite
# Usage: generate_satellite_testcase TLE_file QTH_file satellite_name start_time track_time output_filename
# \param TLE_file File containing TLE data (should not contain more than approx. 20 satellites due to internal restrictions in predict)
# \param QTH_file File containing the QTH
# \param satellite_name Name of satellite that is to be tracked (will use corresponding satellite in the TLE file)
# \param start_time Start time in format e.g. 2015-09-21 11:00
# \param track_time Total tracked time in number of seconds from the start time
function generate_satellite_testcase(){
	tle_file="$1"
	qth_file="$2"
	db_file="$3"
	satellite_name="$4"
	start_time="$5"
	tot_secs="$6"
	testcase_filename="../data/sat_${satellite_name}_$(echo $start_time | sed -r 's/[-|_|:| ]//g').test"

	#move db file where predict can find it
	mkdir -p .predict
	cp $db_file .predict/predict.db
	prevhome="$HOME"
	HOME="."

	#parse tle and qth information into testcase file
	echo "[tle]" > $testcase_filename
	grep -A 2 "$satellite_name" $tle_file | tail -2 >> $testcase_filename

	echo "" >> $testcase_filename
	echo "$(get_qth_string $qth_file)" >> $testcase_filename

	echo "freq=0" >> $testcase_filename

	#parse alat, alon into testcase file
	alatalon=$(echo $(grep "$satellite_name" "$db_file" -A 2 | sed '3!d' | sed -rn 's/(.*)/\1/p'))
	if [ -z "$alatalon" ]; then
		alatalon="No alat, alon" #in case there was no entry for the satellite in the .db file
	fi
	alatalon=$(echo $alatalon | sed -r 's/, alon/, No alon/g')
	echo "alat=$(echo $alatalon | sed -r 's/(.*), .*/\1/g')" >> $testcase_filename
	echo "alon=$(echo $alatalon | sed -r 's/.*, (.*)/\1/g')" >> $testcase_filename

	#predict orbit
	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	faketime "$start_time" predict -t $tle_file -q $qth_file -s &
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_SAT $satellite_name"))
		echo ${predict_response[@]} >> "test_file.dat"
		time=$(./predict_client "GET_TIME")
		doppler_shift=($(./predict_client "GET_DOPPLER $satellite_name"))
		satname=${predict_response[0]}
		lon=$(echo "360-${predict_response[1]}" | bc) #convert to N/E
		lat=${predict_response[2]}
		az=${predict_response[3]}
		el=${predict_response[4]}
		aostime=${predict_response[5]}
		footprint=${predict_response[6]}
		range=${predict_response[7]}
		alt=${predict_response[8]}
		vel=${predict_response[9]}
		revolutions=${predict_response[10]}
		visibility=${predict_response[11]}

		case $visibility in
			N)
				visibility="0"
				;;
			D)
				visibility="1"
				;;
			V)
				visibility="2"
				;;
		esac

		phase=${predict_response[12]}
		eclipse_depth=${predict_response[13]}
		squint=${predict_response[14]}
		echo "$time,$lat,$lon,$alt,$az,$el,$doppler_shift,$squint,$phase,$revolutions,$footprint,$range,$vel,$visibility,$eclipse_depth" >> $testcase_filename
	done
	killall predict
	sleep 1
	HOME="$prevhome"
}

# Generate testcase data for sun tracking. 
# Usage: generate_sun_testcase start_time track_time QTH_file output_filename
# \param start_time Start time in format e.g. 2015-09-21 11:00
# \param track_time Total tracked time in number of seconds from the start time
# \param QTH_file File containing the QTH
function generate_sun_testcase(){
	start_time="$1"
	tot_secs="$2"
	qth_file="$3"
	testcase_filename="../data/sun_$(echo $start_time | sed -r 's/[-|_|:| ]//g').test"

	#parse qth information into testcase file
	echo "$(get_qth_string $qth_file)" >> $testcase_filename

	faketime "$start_time" predict -q $qth_file -s &

	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_SUN"))
		time=$(./predict_client "GET_TIME")
		az=${predict_response[0]}
		el=${predict_response[1]}
		dec=${predict_response[2]}
		gha=${predict_response[3]}
		ra=${predict_response[4]}
		echo $time,$az,$el,$dec,$gha,$ra >> $testcase_filename
	done
	killall predict
	sleep 1
}

# Generate testcase data for moon tracking. 
# Usage: generate_moon_testcase start_time track_time QTH_file output_filename
# \param start_time Start time in format e.g. 2015-09-21 11:00
# \param track_time Total tracked time in number of seconds from the start time
# \param QTH_file File containing the QTH
function generate_moon_testcase(){
	start_time="$1"
	tot_secs="$2"
	qth_file="$3"
	testcase_filename="../data/moon_$(echo $start_time | sed -r 's/[-|_|:| ]//g').test"

	#parse qth information into testcase file
	echo "$(get_qth_string $qth_file)" >> $testcase_filename

	faketime "$start_time" predict -q $qth_file -s &

	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_MOON"))
		time=$(./predict_client "GET_TIME")
		az=${predict_response[0]}
		el=${predict_response[1]}
		dec=${predict_response[2]}
		gha=${predict_response[3]}
		ra=${predict_response[4]}
		echo $time,$az,$el,$dec,$gha,$ra >> $testcase_filename
	done
	killall predict
	sleep 1
}

killall predict

#compile and prepare UDP client
gcc -o predict_client predict_client.c

#Elliptic orbit
#Is SDP4, has defined alon, alat in .db-file for squint angle calc
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "MOLNIYA_1-29" "2015-09-26 18:00" "20"

#Geostationary
#Is SDP4, has defined alon, alat in .db-file for squint angle calc
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "THOR_III" "2015-09-26 18:00" "20"

#Sun-synchronous
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "HINODE" "2015-09-26 18:00" "20"

#Frozen orbit
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "ERS-1" "2015-09-26 18:00" "20"

#High Earth Orbit
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "VELA-1" "2015-09-26 18:00" "20"

#Medium Earth Orbit
#Is SDP4, has defined alon, alat in .db-file for squint angle calc
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "GPS_BIIA-10" "2015-09-26 18:00" "20"

#Tundra orbit
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "SIRIUS-1" "2015-09-26 18:00" "20"

#Low Earth Orbit
generate_satellite_testcase "testcase.tle" "testcase.qth" "testcase.db" "ISS" "2015-09-26 18:00" "20"

#sun and moon
generate_sun_testcase "2015-09-20 19:33" "20" "testcase.qth"
generate_sun_testcase "2015-09-21 06:00" "20" "testcase.qth"
generate_moon_testcase "2015-09-20 10:00" "20" "testcase.qth"
generate_moon_testcase "2015-09-20 16:00" "20" "testcase.qth"
