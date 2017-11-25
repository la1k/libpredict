#!/bin/bash

#Obtain TLE-s containing both geosynchronous and non-geosynchronous objects,
#and a list over satellite numbers that should correspond to geosynchronous objects.

#Script was tailored specifically for obtaining list over geostationary satellite
#numbers and TLEs at 2017-04-29. Will need tweaking when running this script far into the future,
#as URLs might change and the PDF containing recent orbit classification has to be re-obtained from
#some source.

#celestrak TLEs
wget https://www.celestrak.com/NORAD/elements/tle-new.txt
wget https://www.celestrak.com/NORAD/elements/weather.txt
wget https://www.celestrak.com/NORAD/elements/noaa.txt
wget https://www.celestrak.com/NORAD/elements/goes.txt
wget https://www.celestrak.com/NORAD/elements/resource.txt
wget https://www.celestrak.com/NORAD/elements/sarsat.txt
wget https://www.celestrak.com/NORAD/elements/dmc.txt
wget https://www.celestrak.com/NORAD/elements/tdrss.txt
wget https://www.celestrak.com/NORAD/elements/argos.txt
wget https://www.celestrak.com/NORAD/elements/geo.txt
wget https://www.celestrak.com/NORAD/elements/intelsat.txt
wget https://www.celestrak.com/NORAD/elements/ses.txt
wget https://www.celestrak.com/NORAD/elements/iridium-NEXT.txt
wget https://www.celestrak.com/NORAD/elements/iridium.txt
wget https://www.celestrak.com/NORAD/elements/orbcomm.txt
wget https://www.celestrak.com/NORAD/elements/globalstar.txt
wget https://www.celestrak.com/NORAD/elements/amateur.txt
wget https://www.celestrak.com/NORAD/elements/x-comm.txt
wget https://www.celestrak.com/NORAD/elements/other-comm.txt
wget https://www.celestrak.com/NORAD/elements/molniya.txt
wget https://www.celestrak.com/NORAD/elements/raduga.txt
wget https://www.celestrak.com/NORAD/elements/gorizont.txt
wget https://www.celestrak.com/NORAD/elements/gps-ops.txt
wget https://www.celestrak.com/NORAD/elements/glo-ops.txt
wget https://www.celestrak.com/NORAD/elements/beidou.txt
wget https://www.celestrak.com/NORAD/elements/galileo.txt
wget https://www.celestrak.com/NORAD/elements/sbas.txt
wget https://www.celestrak.com/NORAD/elements/nnss.txt
wget https://www.celestrak.com/NORAD/elements/musson.txt
wget https://www.celestrak.com/NORAD/elements/science.txt
wget https://www.celestrak.com/NORAD/elements/geodetic.txt
wget https://www.celestrak.com/NORAD/elements/education.txt
wget https://www.celestrak.com/NORAD/elements/engineering.txt
wget https://www.celestrak.com/NORAD/elements/radar.txt
wget https://www.celestrak.com/NORAD/elements/military.txt
wget https://www.celestrak.com/NORAD/elements/cubesat.txt
wget https://www.celestrak.com/NORAD/elements/other.txt

mkdir -p tles
mv *.txt tles/

#geostationary satellites, as defined by celestrak.com
cp tles/geo.txt ../data/geostationary.tle

#all acquired TLEs from above
cat tles/*.txt > ../data/large-tle-collection.tle

#Obtain list over geostationary satellites, approximately up to date at time of acquisition. Will probably not
#contain geosynchronous satellites in general, or inactive satellites in geosynchronous orbits.
#URL might be exact in the future and the format might change, this is kept here for future reference.
wget http://www.satellite-calculations.com/Satellite/satellitelist.php -O satlist.html
html2text satlist.html | sed -r 's/\|/ /g' | awk '{print $3}' | sed -rn 's/^([0-9]+)$/\1/p' > ../data/geosynchronous_satellite_numbers.dat

#obtain list of geosynchronous objects from ESA geosynchronous classifiction report ("Classification of geosynchronous objects", ESA. Issued regularly)
#List is up to date at 2017-01-01 (issue 19). Will not contain later launches or orbital changes.
#PDF file is not included in the repository, will have to be downloaded manually.
pdftotext -f 42 -l 150 ESA_GEO_Classification_Report_issue_19.pdf - | sed -rn 's/^([0-9]+)$/\1/p' >> ../data/geosynchronous_satellite_numbers.dat
