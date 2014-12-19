function HDL=whdr;
% whdr.m
%          bulds 15 header in MODS file in an interactive mode.
%
%

%  12/5/2003           Edited by Hyun-Sook Kim
%--------------------------------------------------------------------

%--- start to work
disp('      ************    Building MODS Headers   ************** ')
disp(' ')
disp(' ')

%--- line 1:                    
disp(' *** 1st header:')   
hdr1=input('    1.  Enter title:    ','s');
    HDL{1} = sprintf(' title = %s',hdr1);

%--- line 2:                    
disp(' *** 2nd header:')   
hdr2=input('    2. Enter # of stations: ');
    HDL{2} = sprintf(' stations = %04.0f',hdr2);
    
%--- line 3:                    
disp(' *** 3rd header:')
hdr3=input('    3. Enter a Start Time: [year,month,day,hour,min,sec]  ');
    jd=julian(hdr3)-2440000;
    dn=datenum(hdr3);
        mon    = datestr(dn,'mmm');
        day    = datestr(dn,'dd');
        year   = datestr(dn,'yyyy');
        hhmmss = datestr(dn,'HH:MM:SS');
     HDL{3} = sprintf(' str_time = %10.4f, %s %s %s %s',jd,mon,day,year,hhmmss);
                    clear jd dn mon day year hhmmss
%--- line 4:                    
disp(' *** 4th header:')
hdr4=input('    4. Enter a End Time: [year,month,day,hour,min,sec]  ');
    jd=julian(hdr4)-2440000;
    dn=datenum(hdr4);
        mon    = datestr(dn,'mmm');
        day    = datestr(dn,'dd');
        year   = datestr(dn,'yyyy');
        hhmmss = datestr(dn,'HH:MM:SS');
     HDL{4} = sprintf(' end_time = %10.4f, %s %s %s %s',jd,mon,day,year,hhmmss);
                    clear jd dn mon day year hhmmss
%--- line 5:                    
disp(' *** 5th header:')
disp('          5. Julian Day Offset - no input required')
    HDL{5} = sprintf(' Jday_offset = 2440000');

%--- line 6:                    
disp(' *** 6th header:')
hdr6=input('    6. Enter the minimum longitude:    ');
     HDL{6} = sprintf(' lng_min = %07.4f',hdr6);

%--- line 7:                    
disp(' *** 7th header:')
hdr7=input('    7. Enter the maximum longitude:    ');
     HDL{7} = sprintf(' lng_max = %07.4f',hdr7);

%--- line 8:                    
disp(' *** 8th header:')
hdr8=input('    8. Enter the minimum latitude:    ');
     HDL{8} = sprintf(' lat_min = %07.4f',hdr8);

%--- line 9:                    
disp(' *** 9th header:')
hdr9=input('    9. Enter the maximum latitude:    ');
     HDL{9} = sprintf(' lat_max = %07.4f',hdr9);
     
%--- line 10:                         
 	HDL{10} = ' format = ascii, record interleaving';
%--- line 11:                         
	HDL{11} = ' type = CTD';
%--- line 12:                         
	HDL{12} = ' fields = depth, temperature, salinity';
%--- line 13:                         
	HDL{13} = ' units = meter, Celcius, PSU';
%--- line 14:                         
  	HDL{14} = sprintf(' creation_date = %s %s %s %s %s',datestr(now,8),datestr(now,3),datestr(now,7),datestr(now,13),datestr(now,10));
%--- line 15:                         
    HDL{15} = 'END';

%--------------------------------------------------------------------
% end of writing headers
%--------------------------
disp(' ')
disp(' ')
disp('          ----- Finish Setting Up Headers of MODS File -----  ')

