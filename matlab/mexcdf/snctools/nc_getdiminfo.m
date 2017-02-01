function dinfo = nc_getdiminfo ( arg1, arg2 )
% NC_GETDIMINFO:  returns metadata about a specific NetCDF dimension
%
% DINFO = NC_GETDIMINFO(NCFILE,DIMNAME) returns information about the
% dimension DIMNAME in the netCDF file NCFILE.
%
% DINFO = NC_GETDIMINFO(NCID,DIMID) returns information about the
% dimension with numeric ID DIMID in the already-opened netCDF file
% with file ID NCID.  This form is not recommended for use from the
% command line.
%
% Upon output, DINFO will have the following fields.
%
%    Name:  
%        a string containing the name of the dimension.
%    Length:  
%        a scalar equal to the length of the dimension
%    Unlimited:  
%        A flag, either 1 if the dimension is an unlimited dimension
%        or 0 if not.
%
% In case of an error, an exception is thrown.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_getdiminfo.m 2659 2009-04-01 17:38:36Z johnevans007 $
% $LastChangedDate: 2009-04-01 13:38:36 -0400 (Wed, 01 Apr 2009) $
% $LastChangedRevision: 2659 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error(nargchk(2,2,nargin,'struct'));
error(nargoutchk(1,1,nargout,'struct'));

backend = snc_read_backend(arg1);
switch(backend)
	case 'tmw'
		dinfo = nc_getdiminfo_tmw(arg1,arg2);
	case 'java'
		dinfo = nc_getdiminfo_java(arg1,arg2);
	case 'mexnc'
		dinfo = nc_getdiminfo_mexnc(arg1,arg2);
	otherwise
		error('SNCTOOLS:nc_getdiminfo:unhandledBackend', ...
		      'Unhandled backend ''%s''', backend );
end













