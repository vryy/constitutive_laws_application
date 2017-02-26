#!/bin/sh
gfortran -shared -fPIC -o libPlaxis.so UdsmDP.for usrlib.for

