function status = whydro (hydnam,header,hinfo,htype,z,t,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%  function status = whydro (hydnam,header,hinfo,htype,z,t,s);
%
%  Write hydrography data in MODS ASCII format.                      %
%
%  On Input:
%
%     hydnam  input hydrography data filename (string).
%     header  file header text matrix.
%     hinfo   Cast header information matrix
%             hinfo(nsta,1):     nhvar   - number of variables
%             hinfo(nsta,2):     nhpts   - number of data points
%             hinfo(nsta,3):     castid  - cast identifier
%             hinfo(nsta,4):     hlng    - longitude
%             hinfo(nsta,5):     hlat    - latitude
%             hinfo(nsta,6):     hdpth   - maximum depth
%             hinfo(nsta,7):     htime   - time
%             hinfo(nsta,.):     hscl(1:nhvar) - data scales
%             hinfo(nsta,..):    hflag   - special flag
%    htype    hydrography type identification.
%    z        first  hydrographic field, usually depth data.
%    t        second hydrographic field, usually temperature data.
%    s        third  hydrographic field, usually salinity data.
%
%  On Output:
%
%    status   Exit status flag.  [0] no errors.  [other] error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------------
%  Open hydrography file.
%-------------------------------------------------------------------------------

[hid,message] = fopen (hydnam,'w');

if (hid<0),
   disp (' ');
   disp ('***Error:  WHYDRO - Cannot open file:');
   disp (['           ',setstr(34),hydnam,setstr(34)]);
   disp (['           ',setstr(34),message,setstr(34)]);
   disp (' ');
   status = hid;
   return;
end;

%-------------------------------------------------------------------------------
%  Write header.  Extract number of stations.
%-------------------------------------------------------------------------------

[nline nchar] = size (header);

for n = 1:nline,

   ind = find (abs(header(n,:))~=0);

   count = fprintf (hid,'%s\n',header(n,ind));

   [message status] = ferror (hid);
   if (status~=0),
      disp (' ');
      disp (['***Error:  WHYDRO - Error writing header line ',num2str(n), ...
                                                          ' to file:']);
      disp (['           ',setstr(34),hydnam,setstr(34)]);
      disp (['           ',setstr(34),message,setstr(34)]);
      disp (' ');
      return;
   end;

end;

%-------------------------------------------------------------------------------
% Read hydrographic data.
%-------------------------------------------------------------------------------

[nsta npts] = size (z);

for n=1:nsta

   nhvar = hinfo(n,1);
   nhpts = hinfo(n,2);

   if (nhvar>2),
      status = write_cast (hid,hinfo(n,1:(8+nhvar))',htype(n,:), ...
                                  z(n,1:nhpts),t(n,1:nhpts),s(n,1:nhpts));
     else
      status = write_cast (hid,hinfo(n,1:(8+nhvar))',htype(n,:), ...
                                  z(n,1:nhpts),t(n,1:nhpts),t(n,1:nhpts));
   end;

   if(status~=0)
      disp (' ');
      disp (['***Error:  WHYDRO - unable to write cast number ',num2str(n), ...
                                ' (id=',num2str(hinfo(n,3)),') to file:']);
      disp (['           ',setstr(34),hydnam,setstr(34)]);
      disp (' ');
      return;
   end;

end;

%-------------------------------------------------------------------------------
% Close hydrography file.
%-------------------------------------------------------------------------------

status = fclose (hid);

if (status<0),
   disp (' ');
   disp ('***Error:  WHYDRO - Cannot close file:');
   disp (['           ',setstr(34),hydnam,setstr(34)]);
   disp (' ');
   return;
end;
