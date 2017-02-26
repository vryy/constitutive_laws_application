#!/bin/sh
reset
gfortran -shared -fPIC -o libUmat.so CACAMA.f sdvini.f

