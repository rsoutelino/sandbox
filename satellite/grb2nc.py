#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
import glob
################################################################################

filelist = glob.glob("*.grb")
filelist.sort()

for File in filelist:
    print "Converting %s\n" %(File) 
    os.system("cdo -f nc copy %s %s.nc" %(File, File[:-4]))
    os.system("rm %s" %File)


