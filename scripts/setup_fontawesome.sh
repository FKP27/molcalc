#!/bin/bash

cd molcalc/static

VERSION=5.15.4

wget https://github.com/FortAwesome/Font-Awesome/releases/download/$VERSION/fontawesome-free-$VERSION-web.zip

unzip fontawesome-free-$VERSION-web.zip

rm -r fontawesome
mv fontawesome-free-$VERSION-web fontawesome
rm fontawesome-free-$VERSION-web.zip
