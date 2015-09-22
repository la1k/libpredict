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
	satellite_name="$3"
	start_time="$4"
	tot_secs="$5"
	testcase_filename="testcase_satellite_${satellite_name}_$(echo $start_time | sed -r 's/[-|_|:| ]//g').dat"

	#parse tle and qth information into testcase file
	echo "[tle]" > $testcase_filename
	grep -A 2 $satellite_name $tle_file | tail -2 >> $testcase_filename

	echo "" >> $testcase_filename
	echo "$(get_qth_string $qth_file)" >> $testcase_filename

	echo "freq=0" >> $testcase_filename

	#predict orbit
	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	faketime "$start_time" predict -t $tle_file -q $qth_file -s &
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_SAT $satellite_name"))
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
		orbitnum=${predict_response[10]}
		visibility=${predict_response[11]}
		phase=${predict_response[12]}
		eclipse_depth=${predict_response[13]}
		squint=${predict_response[14]}
		echo "$time,$lat,$lon,$alt,$az,$el,$doppler_shift" >> $testcase_filename
	done
	killall predict
	sleep 1
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
	testcase_filename="testcase_sun_$(echo $start_time | sed -r 's/[-|_|:| ]//g').dat"

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
		echo $time,$az,$el >> $testcase_filename
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
	testcase_filename="testcase_moon_$(echo $start_time | sed -r 's/[-|_|:| ]//g').dat"

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
		echo $time,$az,$el >> $testcase_filename
	done
	killall predict
	sleep 1
}

killall predict

#compile and prepare UDP client
gcc -o predict_client predict_client.c

#satellites
generate_satellite_testcase "testcase.tle" "testcase.qth" "OSCAR-7" "2015-09-20 19:15" "20"
generate_satellite_testcase "testcase.tle" "testcase.qth" "OSCAR-7" "2015-09-20 19:31" "20"

#sun and moon
generate_sun_testcase "2015-09-20 19:33" "20" "testcase.qth"
generate_sun_testcase "2015-09-21 06:00" "20" "testcase.qth"
generate_moon_testcase "2015-09-20 10:00" "20" "testcase.qth"
generate_moon_testcase "2015-09-20 16:00" "20" "testcase.qth"
