clear all;close all;clc

for k = 33:92
  eval(['!mv QS_STGRD3_20050',num2str(k),'* QS_STGRD3_20050',num2str(k)]);
end

