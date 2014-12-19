function [header,hinfo,htype,z,t,s]=rhydro(hydnam,fillnan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%  function [header,hinfo,htype,z,t,s]=rhydro(hydnam,fillnan);
%
%  Read and decode file header for hydrography data written in MODS  %
%  ASCII format.                                                     %
%
%  On Input:
%
%     hydnam  input hydrography data filename (string).
%     fillnan   flag for filling empth z,t,s locations with NaN's.
%                  fillnan=0  ->  do NOT fill.       [Default]
%                  fillnan=1  ->  fill with NaN's.
%
%  On Output:
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Set NaN fill test.

   if (nargin<2), fillnan = 0; end;
   testnan = rem(round(fillnan),2);

%  Open hydrography file.

   [hid,message]=fopen(hydnam,'r');
   if (hid<0),
      error(['Cannot open ',setstr(34),hydnam,setstr(34), ...
             ',  ',setstr(34),message,setstr(34),'.']);
   end;

%  Read in header.

   nhdr=1;
   txt=fgetl(hid);
   while (txt(1:3)~='END')
     lenstr=length(txt);
     header(nhdr,1:lenstr)=txt;
     ind_nblk = find (txt~=' ');
     ind_eq = find (txt=='=');
     if ( (ind_eq(1)-ind_nblk(1)+1) == 10),
       if (txt(ind_nblk(1):ind_eq(1))=='stations ='),
         nsta=str2num(txt((ind_eq(1)+1):lenstr));
       end;
     end;
     nhdr=nhdr+1;
     txt=fgetl(hid);
   end
   lenstr=length(txt);
   header(nhdr,1:lenstr)=txt;

% Read hydrographic data.

  for n=1:nsta
     [line,type,d,f1,f2,err,message]=read_cast(hid);
     if(err~=0)
        disp (' ');
        disp ('***Error:  RHYDRO - unable to read cast in file:');
        disp (['           ',setstr(34),hydnam,setstr(34)]);
        if (length(line)>2),
           disp (['           cast number ',num2str(n),'   castid ', ...
                                                num2str(line(3))]);
         else
           disp (['           cast number ',num2str(n)]);
        end;
        disp (['           ',setstr(34),message,setstr(34)]);
        disp (' ');
        error;
     end
     nhvar=line(1);
     nhpts=line(2);
     hinfo(n,1:8+nhvar)=line(1:8+nhvar)';
     lenstr=length(type);
     htype(n,1:lenstr)=type(1:lenstr);
     z(n,1:nhpts)=d(1:nhpts)';
     t(n,1:nhpts)=f1(1:nhpts)';
     if(nhvar==3), s(n,1:nhpts)=f2(1:nhpts)'; end
  end

% Fill empty spots with NaN's,  if so requested.

  if (testnan),
    nhpmx = max(hinfo(:,2));
    for n=1:nsta
      nhpts=hinfo(n,2);
      nhbar=hinfo(n,1);
      if (nhpts<nhpmx),
        z(n,(nhpts+1):nhpmx) = NaN.*ones(size((nhpts+1):nhpmx));
        t(n,(nhpts+1):nhpmx) = NaN.*ones(size((nhpts+1):nhpmx));
        if(nhvar==3),
          s(n,(nhpts+1):nhpmx) = NaN.*ones(size((nhpts+1):nhpmx));
        end;
      end;
    end;
  end;

  return
