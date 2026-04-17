
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
    [age_14c, id2]=read_14c_calibrated(cores{i,1});
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
    depint=diff(depth_14c);
    dur=diff(Bage_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2) depint(ind2)];
    mindur(i)=min(dur(ind2));
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
end

for j = 1:length(cores2)
    [radio_depthage, Bchron_depthage]=func_combine_files(j);
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
    dur_depth=diff(depth_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2) dur_depth(ind2)];
    mindur(i)=min(dur(ind2));
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
end

%asr=all_sed_ratios;
asr=[];
ind=find(cores_meansr>=8 & mindur<=4);
length(ind)
%cores_ind=find(cores_meansr>=2.5 & mindur<6);
%temp_ind=find(cores_meansr>=8 & mindur>=4);
%cores{temp_ind}
%cores_ages{ind}
for i=1:length(ind)
    asr=[asr; cores_dur{ind(i)}; [NaN NaN NaN]];
end

asrnan=asr;

figure(2)

MIN_SPACING=.5; MAX_SPACING=4;
ind=find(asr(:,2)>=MIN_SPACING & asr(:,2)<=MAX_SPACING);
asr=asr(ind,:);

ind=find((asrnan(:,2)>=MIN_SPACING & asrnan(:,2)<=MAX_SPACING) | isnan(asrnan(:,1)) );
asrnan=asrnan(ind,:);


disp('----------------------------------------------------')
disp(['MIN SPACING = ' num2str(MIN_SPACING)])



ratio_bins=[0.15:.1:1];
ratio_bins=[ratio_bins(1:end) fliplr(1./ratio_bins(1:end))];
bin_ctr=geomean([ratio_bins(1:end-1); ratio_bins(2:end)]);
bin_ctr=round(100*bin_ctr)/100;
[n, binID]=histc(asr(:,1),ratio_bins);

% ratio_bins=10.^([-.9:.2:-.1 .1:.2:.9]);
% bin_ctr=geomean([ratio_bins(1:end-1); ratio_bins(2:end)]);
% [n, edges, binID]=histcounts(asr(:,1),ratio_bins);


% This histogram does NOT account for the fact that different sed rate
% observations have different durations
% subplot(211)
% bar(n(1:end-1))
% set(gca,'XTick',[1:16],'XTickLabel',num2str([bin_ctr']))
% axis tight

% Use binID and asr(:,2) to sum durations for each binned ratio
for i=1:length(binID)
    ind=find(binID==i);
    freq(i)=sum(asr(ind,2)); %total number of kyr with this ratio
end
%sum(freq)
% subplot(211)
% %hist(asr(:,2))
% bar(freq(1:end-1))
% set(gca,'XTick',[1:16],'XTickLabel',num2str([bin_ctr']))
% axis tight

duration=round(10*asr(:,2));
dep_int=round(asr(:,3));
%duration(duration==0)=1;
srw=[]; % Each sed rate sample represents 0.1 kyr
srdw=[]; % Each sed rate sample represents 1 cm
srdw_nan=[];
for i=1:length(asr)
    if duration(i)>0
        srw=[srw;asr(i,1)*ones(duration(i),1)];
    end
    if dep_int(i)>0
        srdw=[srdw;asr(i,1)*ones(dep_int(i),1)];
        %dep_int(i)
    end
end
for i=1:length(asrnan)
    if asrnan(i,3)>0
        srdw_nan=[srdw_nan;asrnan(i,1)*ones(round(asrnan(i,3)),1)];
        %round(asrnan(i,3))
    else
        srdw_nan=[srdw_nan; nan];
    end
end
Nsrw=length(srw);
Nsrdw=length(srdw);
disp(['Total number of kyr: ' num2str(Nsrw/10)])
disp(['Total number of cm: ' num2str(Nsrdw)])
mean_srw=mean(srw)
mean_srdw=mean(srdw)

subplot(211)    
[n, edges, binID]=histcounts(srw,ratio_bins);
%n=histc(srw,ratio_bins);
bar(n(1:end))
set(gca,'XTick',[1:length(bin_ctr)],'XTickLabel',num2str([bin_ctr']))
axis tight
title('Sed Rate Ratio')
ylabel('Count per 100 yr')

subplot(212)    
[nd, edges, binID]=histcounts(srdw,ratio_bins);
%nd=histc(srdw,ratio_bins);
bar(nd(1:end))
set(gca,'XTick',[1:length(bin_ctr)],'XTickLabel',num2str([bin_ctr']))
axis tight
ylabel('Count per 1 cm')
xlabel('Sed Rate Ratio')

%%
figure(3)
histogram(log(srw),[-2.242:.065:2])
set(gca,'XTick',[-1.38 -0.0812 0 0.0816 1.38],'XTickLabel',[.25 .922 1 1.085 4])
axis tight
xlabel('log(nSR)')
ylabel('0.1 kyr increments')
title('Reproducing nSR histogram from Lin et al')