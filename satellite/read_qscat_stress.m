function [lon,lat,data]=read_qscat_stress(fn,typ);
%
% A function to read QuikSCAT wind stress HDF files.
% Dimitris Menemenlis <menemenlis@jpl.nasa.gov>
%  USAGE EXAMPLE
%   fn  = 'QS_STGRD3_1999200.03Jan101328';
%   typ = 'asc_stress_Liu_U';
%   [lon,lat,data] = read_qscat_stress(fn,typ);
%
%  OUTPUT
%   lon : longitude
%   lat : latitude
%   data: 2-D output field of type "typ"
%
%  INPUT
%   fn  : file name
%   typ :
%    asc_stress_Liu_U   Liu&Tang   ascending  zonal wind stress (N/m^2)
%    asc_stress_Liu_V   Liu&Tang   ascending  meridional stress (N/m^2)
%    asc_cd_Liu         Liu&Tang   ascending  drag coefficient
%    asc_stress_Large_U Large&Pond ascending  zonal wind stress (N/m^2)
%    asc_stress_Large_V Large&Pond ascending  meridional stress (N/m^2)
%    asc_cd_Large       Large&Pond ascending  drag coefficient
%    asc_rain_flag      ascending  rain flag
%    asc_time_frac      ascending  time in fraction of a day
%    des_stress_Liu_U   Liu&Tang   descending zonal wind stress (N/m^2)
%    des_stress_Liu_V   Liu&Tang   descending meridional stress (N/m^2)
%    des_cd_Liu         Liu&Tang   descending drag coefficient
%    des_stress_Large_U Large&Pond descending zonal wind stress (N/m^2)
%    des_stress_Large_V Large&Pond descending meridional stress (N/m^2)
%    des_cd_Large       Large&Pond descending drag coefficient
%    des_rain_flag      descending rain flag
%    des_time_frac      descending time in fraction of a day
  
  if exist(fn) ~= 2
    error('file name not recognized')
  end

  [jj1 jj2 ii1 ii2 lat1 lat2 lon1 lon2]=coord_qscat_stress;
  lon=(360/(ii2+1))*((ii1:ii2)+.5);
  lat=(180/(jj2+1))*((jj1:jj2)+.5)-90;

  sd_id=hdfsd('start',fn,'read');  
  attr_idx1=hdfsd('findattr',sd_id,'date_of_average');
  [date_ave,status]=hdfsd('readattr',sd_id,attr_idx1);

  for sds_idx=0:16
    if sds_idx==16
      error('typ not recognized')
    end
    sds_id=hdfsd('select',sd_id,sds_idx);
    [dsname,dsndims,dsdims,dstype,dsatts,stat]=hdfsd('getinfo',sds_id);
    if strcmp(typ,dsname)
      break
    end
  end

  start=[jj1 ii1];
  stride = []; 
  edges = [(jj2-jj1) (ii2-ii1)]+1; 
  [data,status] = hdfsd('readdata',sds_id,start,stride,edges);
  attr_idx=0;
  [scale_factor,status]=hdfsd('readattr',sds_id,0);
  [add_offset,status]=hdfsd('readattr',sds_id,2);
  data=double(data)*double(scale_factor)-double(add_offset);
  status = hdfsd('endaccess',sds_id);
  status = hdfsd('end',sd_id);
