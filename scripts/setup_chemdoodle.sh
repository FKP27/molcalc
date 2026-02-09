#!/bin/bash

cd molcalc/static/

wget https://web.chemdoodle.com/downloads/ChemDoodleWeb-11.0.0.zip

unzip ChemDoodleWeb-11.0.0.zip

mkdir chemdoodleweb

cp -r ChemDoodleWeb-11.0.0/install/* chemdoodleweb/
cp -r ChemDoodleWeb-11.0.0/src/* chemdoodleweb/

rm -r ChemDoodleWeb-11.0.0*

