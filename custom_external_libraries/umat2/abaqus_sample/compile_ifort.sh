#!/bin/sh
ifort -shared -fPIC -o libPlaxis.so UdsmDP.for usrlib.for

