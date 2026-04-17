clear cores
cores={'MD84-527';'MD88-770';'SO42-74KL';'EW9209-1JPC';'EW9209-2JPC';...
    'EW9209-3JPC';'GeoB7920-2';'GeoB9508-5';'GeoB9526';'GIK13289-2';...
    'KF13';'KNR31-GPC5all';'MD03-2698';'MD99-2334';'SU81-18';'H214';...
    'MD01-2421';'ODP1145';'RS147-07';'SO50-31KL';...
    'TR163-22';'V19-30';'W8709A-8';'MD02-2589';'PS2489-2comp';...
    'RC16-84';'V24-253';'M35003-4';'MD99-2339';'OCE205-100GGC';...
    'PO200-10-6-2';'AHF16832';'DSDP594';'EW9504-03';'EW9504-04';...
    'EW9504-05';'EW9504-08';'EW9504-09';'GIK17961-2';'GIK17964-2';...
    'MD97-2120';'MD97-2151';'GeoB1711';'GeoB1720-2';...
    'GeoB3202';'KNR159-36'};  
% KH94-3-8, KT90-9-5, V19-28, KT90-9-21 %excluded b/c only 2 14c dates

cores2={'MD95-2042';'MD01-2416';'MD02-2489';'MD98-2181';'W8709A-13';...
    'MD07-3076';'RC11-83';'GIK17940-2';'V35-5'};

all_sed_ratios=[];

for i = 1:length(cores)
    % read in 14C depths and calibrated ages for each core
    [depth_14c, id1]=read_radiocarbon(cores{i,1});
    [age_14c_2p5, age_14c, age_14c_97p5, id2]=read_14c_calibrated(cores{i,1});
    cores_14c{i}=[depth_14c age_14c];
    [Bdepth,Bage]=read_14c_ranges(cores{i,1});
    
    % Bchron age models frequently don't include first and last 14c dates
    % due to even 1-kyr increments of output
    % Add these ages back in
    if depth_14c(1) < min(Bdepth)
        Bdepth=[depth_14c(1); Bdepth];
        Bage=[age_14c(1); Bage];
    end
    if depth_14c(end) > max(Bdepth)
        Bdepth=[Bdepth; depth_14c(end)];
        Bage=[Bage; age_14c(end)];
    end
    
    %Bchron model ages for depth of each 14c measurement
    Bage_14c=interp1(Bdepth,Bage,depth_14c);
    cores_Bage{i}=Bage_14c;
    cores_depth{i}=depth_14c;
    
    %Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
    cores_sr{i}=diff(depth_14c)./diff(Bage_14c);
    cores_ages{i}=[min(Bage_14c) max(Bage_14c)];
    
	% mean sed rate for core
    cores_meansr(i)=(max(depth_14c)-min(depth_14c))/(max(Bage_14c)-min(Bage_14c));
    
    % Divide by mean sed rate to convert to sed rate ratio
	sed_rate_ratio=cores_sr{i}/cores_meansr(i);
    ind1=find(~isnan(sed_rate_ratio));
    
    % Duration of sed rate
    dur=diff(Bage_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2)];
    mindur(i)=min(dur(ind2));
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
end

for j = 1:length(cores2)
    [radio_depthage, Bchron_depthage]=func_combine_files_copy(j);
    i=length(cores)+1;
    cores{i,1}=cores2{j};
    cores_14c{i}=radio_depthage;
    depth_14c=radio_depthage(:,1);
    age_14c=radio_depthage(:,2);
    Bdepth=Bchron_depthage(:,1);
    Bage=Bchron_depthage(:,2);
    
    % Bchron age models frequently don't include first and last 14c dates
    % due to even 1-kyr increments of output
    % Add these ages back in
    if depth_14c(1) < min(Bdepth)
        Bdepth=[depth_14c(1); Bdepth];
        Bage=[age_14c(1); Bage];
    end
    if depth_14c(end) > max(Bdepth)
        Bdepth=[Bdepth; depth_14c(end)];
        Bage=[Bage; age_14c(end)];
    end
    
    %Bchron model ages for depth of each 14c measurement
    Bage_14c=interp1(Bdepth,Bage,depth_14c);
    cores_Bage{i}=Bage_14c;
    cores_depth{i}=depth_14c;

    %Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
    cores_sr{i}=diff(depth_14c)./diff(Bage_14c);
    
	% mean sed rate for core
    cores_meansr(i)=(max(depth_14c)-min(depth_14c))/(max(Bage_14c)-min(Bage_14c));
    cores_ages{i}=[min(Bage_14c) max(Bage_14c)];

    % Divide by mean sed rate to convert to sed rate ratio
	sed_rate_ratio=cores_sr{i}/cores_meansr(i);
    ind1=find(~isnan(sed_rate_ratio));
    
    % Duration of sed rate
    dur=diff(Bage_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2)];
    mindur(i)=min(dur(ind2));
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
end

%asr=all_sed_ratios;
asr=[];
indcores=find(cores_meansr>=8 & mindur<=4);
%cores_ages{ind}
for i=1:length(indcores)
    asr=[asr; cores_dur{indcores(i)}];
end
figure(2)

MIN_SPACING=0.5; MAX_SPACING=4;
ind=find(asr(:,2)>=MIN_SPACING & asr(:,2)<=MAX_SPACING);
asr=asr(ind,:);
disp('----------------------------------------------------')
disp(['MIN SPACING = ' num2str(MIN_SPACING)])
disp(['MAX SPACING = ' num2str(MAX_SPACING)])
disp(['Number of cores = ' num2str(length(indcores))])

ratio_bins=[0.25:.1:1];
ratio_bins=[ratio_bins fliplr(1./ratio_bins)];
bin_ctr=mean([ratio_bins(1:end-1); ratio_bins(2:end)]);
bin_ctr=round(10*bin_ctr)/10;
[n, binID]=histc(asr(:,1),ratio_bins);

% This histogram does NOT account for the fact that different sed rate
% observations have different durations
% subplot(211)
% bar(n(1:end-1))
% set(gca,'XTick',[1:16],'XTickLabel',num2str([bin_ctr']))
% axis tight

% Use binID and asr(:,2) to sum durations for each binned ratio
for i=1:length(binID)
    ind=find(binID==i);
    freq(i)=sum(asr(ind,2));  %total number of kyr with this ratio
end
%sum(freq)
% subplot(211)
% %hist(asr(:,2))
% bar(freq(1:end-1))
% set(gca,'XTick',[1:16],'XTickLabel',num2str([bin_ctr']))
% axis tight

duration=round(10*asr(:,2));
%duration(duration==0)=1;
srw=[]; % Each sed rate sample represents 0.1 kyr
for i=1:length(asr)
    if duration(i)>0
        srw=[srw;asr(i,1)*ones(duration(i),1)];
    end
end
Nsrw=length(srw);
disp(['Total number of kyr: ' num2str(Nsrw/10)])
mean_srw=mean(srw)

%subplot(212)    
n=histc(srw,ratio_bins);
bar(n(1:end-1))
set(gca,'XTick',[1:16],'XTickLabel',num2str([bin_ctr']))
axis tight

srws=sort(srw);
% disp('ratio distribution: (2.5% 50% 97.5%)')
% disp(num2str([srws(round(0.025*Nsrw)) srws(round(0.5*Nsrw)) srws(round(0.975*Nsrw))]))

ctr1=min(find(srws>1));
ctr2=max(find(srws<1));
disp('Percent <1, 1, >1')
disp(num2str(100*[ctr2 (ctr1-ctr2) (Nsrw-ctr1)]./Nsrw))
disp('Percent <0.9, .9-1.1,  >1.1')
disp(num2str(100*[min(find(srws>.9)) min(find(srws>1.1))-min(find(srws>.9))...
    Nsrw-min(find(srws>1.1))]./Nsrw))
disp('Percent <0.82, .82-1.22,  >1.22')
disp(num2str(100*[min(find(srws>.82)) min(find(srws>1.22))-min(find(srws>.82))...
    Nsrw-min(find(srws>1.22))]./Nsrw))
disp('Percent <0.75, .75-1.33,  >1.33')
disp(num2str(100*[min(find(srws>.75)) min(find(srws>1.33))-min(find(srws>.75))...
    Nsrw-min(find(srws>1.33))]./Nsrw))
disp('Percent <0.67, .67-1.5,  >1.5')
disp(num2str(100*[min(find(srws>.67)) min(find(srws>1.5))-min(find(srws>.67))...
    Nsrw-min(find(srws>1.5))]./Nsrw))
disp('Percent <0.5,  >2')
disp(num2str(100*[min(find(srws>.5)) (Nsrw-min(find(srws>2)))]./Nsrw))
disp(' ')
disp('ratio for (1:1 - 34%) and (1:1 + 34%)')
disp(num2str([srws(round(.317*ctr1)) srws(ctr2+round(.683*(Nsrw-ctr2)))]))
% disp('ratio distribution: (15.85% 50% 84.15%)')
% disp(num2str([srws(round(0.1585*Nsrw)) srws(round(0.5*Nsrw)) srws(round(0.8415*Nsrw))]))

disp('Percent samples outside of range')
disp(num2str([length(find(srw<.25)) length(find(srw>4))]./Nsrw))

%srws([1:round(Nsrw/11):Nsrw Nsrw])'
