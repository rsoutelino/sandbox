function [jj1, jj2, ii1, ii2,lat1,lat2,lon1,lon2]=coord_qscat_stress()
%
% function needed by read_qscat_stress
% Dimitris Menemenlis <menemenlis@jpl.nasa.gov>  
% Select latitude and longitudes
  lat1=-90; lat2=90;
  if lat1 < -90|lat2 < -90|lat1 > 90|lat2 > 90
    error('ERROR: Latitudes must be between -90 and 90')
  end
  
  % Make sure that lat2 is greater than lat1
  if lat1 > lat2 
    itmp=lat1;
    lat1=lat2;
    lat2=itmp;
  end
  
  % The last grid point is in cell 719.  Reduce lat 90. to 89.9
  if lat2 == 90 
    lat2=89.9;
  end
  
  lon1=0; lon2=360;
  if lon1 < 0|lon2 < 0 | lon1 > 360 | lon2 > 360
    error('ERROR: Longitudes must be between 0 and 360')
  end
  
  % Make sure that lon2 is greater than lon1
  if lon1 > lon2 
    itmp=lon1;
    lon1=lon2;
    lon2=itmp;
  end
  
  % The last grid point is in cell 1439.  Wrapping is not done here,
  % so 360, must be reduced to 359.9.
  if lon2 == 360 
    lon2=359.9;
  end
  
  % Determine grid points from the latitudes and longitudes
  dx=(360./1440.);
  ii1=fix(lon1/dx);
  ii2=fix(lon2/dx);
  
  dy=(180./720.);
  jj1=fix((lat1+90.)/dy);
  jj2=fix((lat2+90.)/dy);
