#!/usr/bin/python
import os

origin = "/media/SAMSUNG/leandro_acer_bkp/"
# origin = "/home/simone/leandro_acer_backup/teste/"
destination = "/home/simone/leandro_acer_backup/apresentacoes/"
pptlist = []
copylist = []

for root, directories, filenames in os.walk(origin):
#	print root
#	print directories	
#	print filenames
	for name in filenames:
	    if "ppt" in name or "odp" in name or "PPT" in name:
			pptlist.append(name)
			copy = "cp '%s' %s" %( os.path.join(root, name), destination )
			copylist.append(copy)
			os.system(copy)
			print "   COPYING  %s\n\n"  %copy
			
		

