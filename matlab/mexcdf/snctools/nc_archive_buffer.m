function new_data = nc_archive_buffer ( input_buffer, ncfile, record_variable )
% NC_ARCHIVE_BUFFER:  Tacks on new data from simple matlab structure to an unlimited-dimension netcdf file
% 
% This function is deprecated.  Please use nc_addnewrecs.m instead.  All 
% this function really does is call nc_addnewrecs anyway.


warning ( 'SNCTOOLS:NC_ARCHIVE_BUFFER:deprecated', ...
    '%s is deprecated and may be removed in a future version of SNCTOOLS.',...
    upper(mfilename));

new_data = [];

error(nargchk(2,3,nargin,'struct'));

switch nargin
case 2
	new_data = nc_addnewrecs ( ncfile, input_buffer );
case { 3, 4 }
	new_data = nc_addnewrecs ( ncfile, input_buffer, record_variable );
end


return;



