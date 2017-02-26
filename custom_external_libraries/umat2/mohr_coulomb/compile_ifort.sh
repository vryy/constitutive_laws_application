#!/bin/sh
reset
ni ifort -shared -fPIC -o libUmat.so umat.for

