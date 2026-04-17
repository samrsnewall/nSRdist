function [depth,age]=read_14c_ranges(core)
 
 filename=['Output/' core '_radiocarbonranges.txt'];
 fid=fopen(filename,'r');
 C=textscan(fid,'%[^\n]',1);
 D=textscan(fid,'%f %f %f %f');  %D{1}=depth, D{3}=age
 fclose(fid);
 depth=D{1};
 age=D{3};

 
