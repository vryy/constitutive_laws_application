#!/bin/sh
reset
gfortran -shared -fPIC -o libUmat.so umat.for umat_hcea.for

