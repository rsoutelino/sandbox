function status = write_cast (fid,line,htype,z,t,s);
%
%function status = write_cast (fid,line,htype,z,t,s);
%
% This function writes the two lines of cast
% header information and the cast itself.
%
% ------
% Input:
% ------
%
% fid.......File identifier for cast data.
% line.....Cast header information.
%             line(1):          nhvar   - number of variables
%             line(2):          nhpts   - number of data points
%             line(3):          castid  - cast identifier
%             line(4):          hlng    - longitude
%             line(5):          hlat    - latitude
%             line(6):          hdpth   - maximum depth
%             line(7):          htime   - time
%             line(8:7+nhvar):  hscl(n) - data scales
%             line(.):          hflag   - special flag
% htype.....hydrography type identification.
% z.........Depth data.
% t.........Temperature data.
% s.........Salinity data.
%
% -------
% Output:
% -------
%
% status....Exit status.  [0] no error  [other] error
%

%-------------------------------------------------------------------------------
% Extract counter information.
%-------------------------------------------------------------------------------

nhvar = line(1);
nhpts = line(2);

%-------------------------------------------------------------------------------
% Write lines of cast header.  Check for error.
%-------------------------------------------------------------------------------

%-----------------------------------
%--- Print constant header data. ---
%-----------------------------------

count = fprintf (fid,' %1d %5d %5d %9.4f %8.4f %7.1f  %10.4f',line(1:7));

[message status] = ferror (fid);
if (status~=0),
   disp(' ');
   disp('***Error:  WRITE_CAST - Unable to write basic cast header data.');
   disp(['           ',setstr(34),message,setstr(34)]);
   return;
end;

%------------------------------
%--- Print variable scales. ---
%------------------------------

for n = 1:nhvar
   count = fprintf (fid,' %8.2e',line(7+n));
   [message status] = ferror (fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write scaling for variable ', ...
                                           num2str(n)]);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

%------------------------------
%--- Print carriage return. ---
%------------------------------

count = fprintf (fid,'\n');

[message status]=ferror(fid);
if (status~=0),
   disp(' ');
   disp(['***Error:  WRITE_CAST - Unable to write carriage return for ', ...
                                             'first line of cast header.']);
   disp(['           ',setstr(34),message,setstr(34)]);
   return;
end;

%--------------------------------
%--- Print unused flag value. ---
%--------------------------------

count = fprintf (fid,' %2d',line(8+nhvar));
[message status]=ferror(fid);
if (status~=0),
   disp(' ');
   disp(['***Error:  WRITE_CAST - Unable to write flag value to cast header.']);
   disp(['           ',setstr(34),message,setstr(34)]);
   return;
end;

%--------------------------------------
%--- Print type string from ' to '. ---
%--------------------------------------

ind = find (htype==setstr(39));
count = fprintf (fid,' %s\n',htype(ind(1):ind(2)));
[message status]=ferror(fid);
if (status~=0),
   disp(' ');
   disp(['***Error:  WRITE_CAST - Unable to write type to cast header.']);
   disp(['           ',setstr(34),message,setstr(34)]);
   return;
end;

%-------------------------------------------------------------------------------
% Convert data and write casts.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------
%--- Determine number of points in full line & number left over. ---
%-------------------------------------------------------------------

nrem  = rem ( nhpts, 10);
nfull = nhpts - nrem;
if (nrem>1),
   indrem = nfull + (1:(nrem-1));
end;

%--------------------
%--- Print depth. ---
%--------------------

if (nfull>0),
   count = fprintf (fid,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n', ...
                                    round(z(1:nfull)./line(8)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write full lines of depth.']);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

if (nrem>1),
   count  = fprintf (fid,'%5d ',round(z(indrem)./line(8)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write leading remaining ', ...
                                       'values of depth.']);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

if (nrem>0),
   count  = fprintf (fid,'%5d\n',round(z(nhpts)./line(8)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write last remaining ', ...
                                       'values of depth.']);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

%--------------------------
%--- Print temperature. ---
%--------------------------

if (nfull>0),
   count = fprintf (fid,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n', ...
                                    round(t(1:nfull)./line(9)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write full lines of ', ...
                                                    'temperature.']);
     disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

if (nrem>1),
   count  = fprintf (fid,'%5d ',round(t(indrem)./line(9)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write leading remaining ', ...
                                       'values of temperature.']);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

if (nrem>0),
   count  = fprintf (fid,'%5d\n',round(t(nhpts)./line(9)) );
   [message status]=ferror(fid);
   if (status~=0),
      disp(' ');
      disp(['***Error:  WRITE_CAST - Unable to write last remaining ', ...
                                       'values of temperature.']);
      disp(['           ',setstr(34),message,setstr(34)]);
      return;
   end;
end;

%-----------------------
%--- Print salinity. ---
%-----------------------

if (nhvar==3),

   if (nfull>0),
      count = fprintf (fid,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n', ...
                                       round(s(1:nfull)./line(10)) );
      [message status]=ferror(fid);
      if (status~=0),
         disp(' ');
         disp(['***Error:  WRITE_CAST - Unable to write full lines of ', ...
                                                          'salinity.']);
         disp(['           ',setstr(34),message,setstr(34)]);
         return;
      end;
   end;

   if (nrem>1),
      count  = fprintf (fid,'%5d ',round(s(indrem)./line(10)) );
      [message status]=ferror(fid);
      if (status~=0),
         disp(' ');
         disp(['***Error:  WRITE_CAST - Unable to write leading remaining ', ...
                                          'values of salinity.']);
         disp(['           ',setstr(34),message,setstr(34)]);
         return;
      end;
   end;

   if (nrem>0),
      count  = fprintf (fid,'%5d\n',round(s(nhpts)./line(10)) );
      [message status]=ferror(fid);
      if (status~=0),
         disp(' ');
         disp(['***Error:  WRITE_CAST - Unable to write last remaining ', ...
                                          'value of salinity.']);
         disp(['           ',setstr(34),message,setstr(34)]);
         return;
      end;
   end;

end,
