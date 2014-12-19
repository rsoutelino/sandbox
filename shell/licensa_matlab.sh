#!/bin/bash -v

# cd /tmp
# /bin/rm /tmp/license.dat
# wget
# http://matlab:'sdf%$t'@dat.cce.usp.br/software/Matlab/Linux/arquivos/license.dat

for mach in macae prainha arraial leblon camboinhas itamambuca conchas itacoatiara grumari buzios urca cabofrio
  do
  scp ../Desktop/R13/license.dat root@$mach:/usr/local/matlab/etc/license.dat
done 
