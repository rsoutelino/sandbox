function [fid,numdat,strdat]= read_header (hydnam);
%
%  function [fid,numdat,strdat]= read_header (hydnam);
%
%  Open an ASCII data file in MODS format.  Read header data.
%
%  ------
%  Input:
%  ------
%
%     HYDNAM  input hydrography data filename (string).
%
%  -------
%  Output:
%  -------
%
%     FID     File identifier.
%     NUMDAT  Numerical data from header.
%                NUMDAT(1):  NSTA - number of stations
%                NUMDAT(2):  STR_TIME - starting time   (modified Julian day)
%                NUMDAT(3):  END_TIME - ending time     (modified Julian day)
%                NUMDAT(4):  JDAY_OFFSET - Julian offset  (day)
%                NUMDAT(5):  LONMIN - minimum longitude
%                NUMDAT(6):  LONMAX - maximum longitude
%                NUMDAT(7):  LATMIN - minimum latgitude
%                NUMDAT(8):  LATMAX - maximum latgitude
%     STRDAT  String data from header.
%                STRDAT(1,:):  File title
%                STRDAT(2,:):  File format
%                STRDAT(3,:):  Probe types
%                STRDAT(4,:):  Data types
%                STRDAT(5,:):  Data units

%-------------------------------------------------------------------------------
%  Initialize all output.
%-------------------------------------------------------------------------------

fid    = [];
numdat = [];
strdat = 'null';

%-------------------------------------------------------------------------------
%  Open hydrography file.
%-------------------------------------------------------------------------------

[fid,message] = fopen (hydnam,'r');
if (fid<0),
   error(['Cannot open ',setstr(34),hydnam,setstr(34), ...
          ',  ',setstr(34),message,setstr(34),'.']);
end;

%-------------------------------------------------------------------------------
%  Read in header.
%-------------------------------------------------------------------------------

txt = fgetl(fid);

while (txt(1:3)~='END')

   lenstr = length (txt);
   ind = findstr(txt,'=');
   type = upper( txt(2:(ind(1)-2)) );

   if strcmp (type,'TITLE')
      strt = ind+1;
      while strcmp(txt(strt),' ')
         strt=strt+1;
      end;
      wkstr = txt(strt:lenstr);
      strdat(1,1:length(wkstr)) = wkstr;
    elseif strcmp (type,'STATIONS')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(1) = str2num(txt(wkind2));
    elseif strcmp (type,'STR_TIME')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1);
      wkind = findstr(txt,',') - 1;
      numdat(2) = str2num(txt(wkind2:wkind));
    elseif strcmp (type,'END_TIME')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1);
      wkind = findstr(txt,',') - 1;
      numdat(3) = str2num(txt(wkind2:wkind));
    elseif strcmp (type,'JDAY_OFFSET')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(4) = str2num(txt(wkind2));
    elseif strcmp (type,'LNG_MIN')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(5) = str2num(txt(wkind2));
    elseif strcmp (type,'LNG_MAX')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(6) = str2num(txt(wkind2));
    elseif strcmp (type,'LAT_MIN')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(7) = str2num(txt(wkind2));
    elseif strcmp (type,'LAT_MAX')
      strt = ind+1;
      wkind = find(~strcmp(txt(strt:lenstr),' '));
      wkind2 = (min(wkind)+strt-1):lenstr;
      numdat(8) = str2num(txt(wkind2));
    elseif strcmp (type,'FORMAT')
      strt = ind+1;
      while strcmp(txt(strt),' ')
         strt=strt+1;
      end;
      wkstr = txt(strt:lenstr);
      strdat(2,1:length(wkstr)) = wkstr;
    elseif strcmp (type,'TYPE')
      strt = ind+1;
      while strcmp(txt(strt),' ')
         strt=strt+1;
      end;
      wkstr = txt(strt:lenstr);
      strdat(3,1:length(wkstr)) = wkstr;
    elseif strcmp (type,'FIELDS')
      strt = ind+1;
      while strcmp(txt(strt),' ')
         strt=strt+1;
      end;
      wkstr = txt(strt:lenstr);
      strdat(4,1:length(wkstr)) = wkstr;
    elseif strcmp (type,'UNITS')
      strt = ind+1;
      while strcmp(txt(strt),' ')
         strt=strt+1;
      end;
      wkstr = txt(strt:lenstr);
      strdat(5,1:length(wkstr)) = wkstr;
   end;

   txt = fgetl(fid);

end
