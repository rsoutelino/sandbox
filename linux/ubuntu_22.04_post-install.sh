USERNAME=rsoutelino
GROUPNAME=rsoutelino
DEFAULT_DIRS='/data /hot /static /scratch /source /config /archive'

# useradd -u 1005 $USERNAME
# echo $USERNAME:password | chpasswd

# groupadd -g 1005 $GROUPNAME
# usermod -G $GROUPNAME -a $USERNAME

mkdir $DEFAULT_DIRS
chown -R $USERNAME:$GROUPNAME $DEFAULT_DIRS
chmod g+rw $DEFAULT_DIRS

echo 'deb http://archive.canonical.com/ubuntu jammy partner' >> /etc/apt/sources.list
apt update

#apt -y install mlocate tilix ntp rsync wget curl python3-pip ipython3 python3-virtualenv python3-virtualenvwrapper git screen vim htop unzip zsh docker.io build-essential gfortran m4 zsh nco ncview vim-gtk subversion libnetcdf-dev libhdf5-serial-dev libkernlib1-gfortran netcdf-bin hdf5-tools libgsl0-dev libgeos-dev libproj-dev python-dev
#pip install seawater pandas pydap netcdf4 pyephem xarray matplotlib ruamel.yaml dask zarr seaborn docker-compose
chmod 777 /var/run/docker.sock
