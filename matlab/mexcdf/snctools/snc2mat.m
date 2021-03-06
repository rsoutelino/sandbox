function snc2mat ( ncfile, matfile ) %#ok<INUSD>
% SNC2MAT:  saves netcdf file to *.mat format.  
%
% SNC2MAT(NCFILE,MATFILE) will save the netCDF file NCFILE to the mat-file 
% MATFILE.  This function is deprecated and may disappear in a future release
% of SNCTOOLS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: snc2mat.m 2663 2009-04-10 19:06:54Z johnevans007 $
% $LastChangedDate: 2009-04-10 15:06:54 -0400 (Fri, 10 Apr 2009) $
% $LastChangedRevision: 2663 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning ( 'SNCTOOLS:snc2mat:deprecated', ...
    'This function is deprecated and may be removed in a future version of SNCTOOLS.' );




%
% create the MATLAB file
ncdata = nc_getall ( ncfile );

fnames = fieldnames ( ncdata );
save_command = '';
global_atts = [];
for j = 1:length(fnames)
	theVar = fnames{j};
	if ( strcmp(theVar,'global_atts' ) )
		global_atts = ncdata.global_atts;
	else
		command = sprintf ( '%s = ncdata.%s;', theVar, theVar );
		eval(command);
		save_command = sprintf ( '%s''%s'',', save_command, theVar );
	end
end
if ~isempty(global_atts)
	save_command = sprintf ( '%s''global_atts''', save_command );
else
	%
	% This chops off a bad comma that's not needed if no global attributes.
	save_command(end) = '';
end
save_command = sprintf ( 'save ( matfile, %s );', save_command );
try
	eval(save_command);
catch %#ok<CTCH>
	error ( 'SNCTOOLS:snc2mat:badCommand', ...
        'Could not execute save command ''%s''', save_command);
end


	

return
