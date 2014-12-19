function [line,htype,z,t,s,error,message]=read_cast(fid)
%
%function [line,htype,z,t,s,error,message]=read_cast(fid)
%
% This function reads the two lines of cast
% header information and the cast itself.
%
% Input:
%
% fid.......File identifier for cast data.
%
% Output:
%
% error....Error flag
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
% message...Error message.
% s.........Salinity data.
% t.........Temperature data.
% z.........Depth data.
%

%-------------------------------------------------------------------------------
% Read line of cast header.  Check for error.
%-------------------------------------------------------------------------------

[line count]=fscanf(fid,'%d %d %d %f %f %f %f',7);
[message error]=ferror(fid);

if (error==0),
   nhvar = line(1);
   nend = 7 + nhvar;
   [line(8:nend) count]=fscanf(fid,'%f',nhvar);
   [message error]=ferror(fid);
 else
   line = [];
end,

if (error==0),
   nend = 8 + nhvar;
   [line(nend) count]=fscanf(fid,'%d',1);
   [message error]=ferror(fid);
end,

if (error==0),
   htype=fgets(fid);
   [message error]=ferror(fid);
end,


%-------------------------------------------------------------------------------
% If still no errors, casts and convert data.
%-------------------------------------------------------------------------------

if (error==0),

%-------------------------------------------------------------------------------
%    Extract number of data points.  Read in data
%-------------------------------------------------------------------------------

   nhpts = line(2);

   [z count]=fscanf(fid,'%d',nhpts);
   [message error]=ferror(fid);

   if (error==0),
      [t count]=fscanf(fid,'%d',nhpts);
      [message error]=ferror(fid);
      if ( (error~=0) & (nhvar<3) & strcmp(upper(message),'AT END-OF-FILE.') )
         error = 0;
      end;
   end

   if ( (error==0) & (nhvar==3) ),
      [s count]=fscanf(fid,'%d',nhpts);
      [message error]=ferror(fid);
      if ( (error~=0) & strcmp(upper(message),'AT END-OF-FILE.') )
         error = 0;
      end;
    elseif (error==0),
      s = NaN.*ones(size(t));
   end

%-------------------------------------------------------------------------------
%    If data successfully read, convert to real values.
%-------------------------------------------------------------------------------

   if (error==0),
      z = line(8).*z;
      t = line(9).*t;
      if (nhvar==3), s = line(10).*s; end,
   end

end,
