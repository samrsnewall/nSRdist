function [age2p5, age50, age97p5, id]=read_14c_calibrated(core)
 
 filename=['Output/' core '_radiocarbonCalibratedRanges.txt'];
 fid=fopen(filename,'r');
 C=textscan(fid,'%[^\n]',1);
 D=textscan(fid,'%s %f %f %f');  %D{3}=mean 14C age
 fclose(fid);
 age2p5 = D{2};
 age50=D{3};
 age97p5 = D{4};
 id=D{1};

