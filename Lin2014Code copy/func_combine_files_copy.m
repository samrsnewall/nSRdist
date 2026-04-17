function [radio_depthage, Bchron_depthage]=func_combine_files_copy(core)

% Combine files for cores that have divided Bchron runs

switch core
    case 1
        cores={'MD95-2042a';'MD95-2042b';'MD95-2042c';'MD95-2042d'};
    case 2
        cores={'MD01-2416a';'MD01-2416b';'MD01-2416c'};
    case 3
        cores={'MD02-2489a';'MD02-2489b';'MD02-2489c'};
    case 4
        cores={'MD98-2181a';'MD98-2181b';'MD98-2181c';'MD98-2181d'};
    case 5
        cores={'W8709A-13a';'W8709A-13b'};
    case 6
        cores={'MD07-3076a';'MD07-3076b';'MD07-3076c';'MD07-3076d'};
    case 7
        cores={'RC11-83a';'RC11-83b';'RC11-83c'};
    case 8
        cores={'GIK17940-2a';'GIK17940-2b';'GIK17940-2c';'GIK17940-2d'};
    case 9
        cores={'V35-5a';'V35-5b';'V35-5c';'V35-5d'};
end
    
id1=[]; id2=[];
depth_14c=[]; age_14c=[];
Bdepth=[]; Bage=[];

for i = 1:length(cores)
    % read in 14C depths and calibrated ages for each core
    [depth1_14c, id1a]=read_radiocarbon(cores{i,1});
    id1=[id1;id1a];
    depth_14c=[depth_14c; depth1_14c];
    
    [age_14c_2p5, age1_14c, age_14c_97p5, id2a]=read_14c_calibrated(cores{i,1});
    id2=[id2;id2a];
    age_14c=[age_14c; age1_14c age_14c_2p5 age_14c_97p5];
    
    [Bdepth1,Bage1]=read_14c_ranges(cores{i,1});
    Bdepth=[Bdepth; Bdepth1];
    Bage=[Bage; Bage1];
end

[id1_uniq,ind1_file,ind1_uniq] = unique(id1,'stable');
[id2_uniq,ind2_file,ind2_uniq] = unique(id2,'stable');
radio_depthage=[depth_14c(ind1_file) age_14c(ind2_file,:)];
Bchron_depthage=[Bdepth Bage];    
    