function [depth,id]=read_radiocarbon(core)
 
 filename=['Input/' core '_radiocarbon.txt'];
 fid=fopen(filename,'r');
 C=textscan(fid,'%s',8);
 D=textscan(fid,'%s %f %f %f %f %f %f %f');  %D{4}=depth
 fclose(fid);
 depth=D{4};
 id=D{1};
 
