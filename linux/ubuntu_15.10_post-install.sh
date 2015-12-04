############################################################################################
# - script for post-install of Ubuntu 15.10 x64 on MSL workstation
# Rafael Soutelino - rsoutelino@gmail.com
############################################################################################
# 
# WHAT SHOULD HAVE BEEN DONE BEFORE RUNNING THIS SCRIPT ---------------------------------
#
# (*) Install Ubuntu from CD-ROM or USB
# (*) Configure network and proxy
# (*) Add contrib and non-free repos to sources.list
# (*) create msluser and metocean groups and give ownerships to home folder
# (*) make sure git is able to clone without password
# (*)
# (*)
# (*)
# (*)
# (*)
###########################################################################################

# adding repositories
wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add -
sudo sh -c 'echo "deb http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google.list'

apt-get update

# system utilities and general purpouse software

apt-get -y install clipit
apt-get -y install inkscape
apt-get -y install google-chrome-stable
apt-get -y install unrar-free
apt-get -y install openssh-server
apt-get -y install build-essential
apt-get -y install gfortran
apt-get -y install docker.io
apt-get -y install zsh
apt-get -y install nco
apt-get -y install ncview
apt-get -y install ntp
apt-get -y install vim-gtk 
apt-get -y install subversion
apt-get -y install vlc
apt-get -y install git-core git
apt-get -y install gobby
apt-get -y install gcolor2
apt-get -y install mysql-workbench
apt-get -y install texlive-latex-base texlive-latex-recommended texlive-latex-extra latex-beamer
wget --no-check-certificate http://install.ohmyz.sh -O - | sh

# scientific libraries

apt-get -y install libnetcdf-dev
apt-get -y install libhdf5-serial-1.8.4 libhdf5-serial-dev 
apt-get -y install libkernlib1-gfortran
apt-get -y install netcdf-bin
apt-get -y install hdf5-tools
apt-get -y install libhdf5-dev

cd /source/
git clone git@github.com:metocean/roms.git
cd roms/docker/netcdf-fortran-4.2
./configure
make
make check
make install

# Python 

apt-get -y install python-dev
apt-get -y install python-scipy
apt-get -y install python-matplotlib
apt-get -y install python-mpltoolkits.basemap
apt-get -y install python-mpltoolkits.basemap-data
apt-get -y install python-pandas
apt-get -y install ipython
apt-get -y install python-pip
apt-get -y install python-setuptools
apt-get -y install python-netcdf
apt-get -y install python-yaml 
apt-get -y install python-mako
apt-get -y install ruamel.yaml

pip install seawater 
pip install pydap 
pip install wx 
pip install netcdf4 
pip install pyephem 
pip install xray







