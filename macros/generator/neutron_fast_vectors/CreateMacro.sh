#!/bin/bash

printf "/generator/add external external:poisson:\"$1\"\n/run/beamOn 1000\n" > $2
