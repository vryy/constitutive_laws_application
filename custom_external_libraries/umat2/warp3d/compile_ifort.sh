#!/bin/sh
reset
ni ifort -shared -fPIC -o libUmat.so umats.f

