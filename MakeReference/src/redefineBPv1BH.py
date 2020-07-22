#!/usr/bin/env python

## Redefine base positions so that beagle v4 can run. 
## Buhm Han, 2017-01-03

## Arguments:
## 1. Marker File
## 2. Output File

import sys, os

def redefineBP(_markerfile, _outfile):

    occupied={}
    with open(_markerfile) as fin, open(_outfile, 'w') as fout:
        for l in fin:
            [rsid, bp, al1, al2]=l.split()
            bp=int(bp)
            while bp in occupied:
                bp+=1
            occupied[bp]="%s %d %s %s\n"%(rsid, bp, al1, al2)
        for key in sorted(occupied):
            fout.write(occupied[key])

    return _outfile


if __name__ == '__main__':


    [markerfile, outfile] = sys.argv[1:3]

    redefineBP(markerfile, outfile)
