flatpak remote-add --if-not-exists fedora
flatpak remote-add --if-not-exists fedora oci+https://registry.fedoraproject.org
flatpak remote-add --if-not-exists flathub https://flathub.org/repo/flathub.flatpakrepo
flatpak remote-add --if-not-exists flathub-beta https://flathub.org/beta-repo/flathub-beta.flatpakrepo
flatpak remote-add --if-not-exists kdeapps https://distribute.kde.org/kdeapps.flatpakrepo
flatpak remote-add --if-not-exists gnome-nightly https://nightly.gnome.org/gnome-nightly.flatpakrepo
flatpak update --appstream
flatpak update


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


dnf -y install mlocate python3-pip ipython3 python3-virtualenv python3-virtualenvwrapper htop zsh podman-docker make automake gcc gcc-c++ kernel-devel gfortran m4 nco ncview netcdf-cxx4 netcdf-cxx4-devel python3-devel

pip3 install seawater pandas pydap netcdf4 pyephem xarray matplotlib ruamel.yaml dask zarr seaborn docker-compose

