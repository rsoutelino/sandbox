#!/usr/bin/env python
import sys
import os

expt = 'phd16'
eddy = "IE"

for zlev in [0, 50, 100, 400]:
    command = "python select_eddy_%s.py '%s' '%s' '%s' " %(expt, expt, zlev, eddy)
    os.system(command)