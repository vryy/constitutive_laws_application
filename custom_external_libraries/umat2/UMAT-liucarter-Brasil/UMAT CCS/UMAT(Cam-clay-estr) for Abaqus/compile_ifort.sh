#!/bin/sh
reset
ni ifort -shared -fPIC -o libUmat.so CACAMA.f sdvini.f

