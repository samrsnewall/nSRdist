% Calculate sed rate ratio for every 1 kyr increments and find
% autocorrelation

cores={'MD84-527';'MD88-770';'SO42-74KL';'EW9209-1JPC';'EW9209-2JPC';...
    'EW9209-3JPC';'GeoB7920-2';'GeoB9508-5';'GeoB9526';'GIK13289-2';...
    'KF13';'KNR31-GPC5all';'MD03-2698';'MD99-2334';'SU81-18';'H214';...
    'KH94-3-8';'KT90-9-5';'MD01-2421';'ODP1145';'RS147-07';'SO50-31KL';...
    'TR163-22';'V19-28';'V19-30';'W8709A-8';'MD02-2589';'PS2489-2comp';...
    'RC16-84';'V24-253';'M35003-4';'MD99-2339';'OCE205-100GGC';...
    'PO200-10-6-2';'AHF16832';'DSDP594';'EW9504-03';'EW9504-04';...
    'EW9504-05';'EW9504-08';'EW9504-09';'GIK17961-2';'GIK17964-2';...
    'KT90-9-21';'MD97-2120';'MD97-2151';'GeoB1711';'GeoB1720-2';...
    'GeoB3202';'KNR159-36'};  %size=50

cores2={'MD95-2042';'MD01-2416';'MD02-2489';'MD98-2181';'W8709A-13';...
    'MD07-3076';'RC11-83';'GIK17940-2';'V35-5'};

all_sed_ratios=[];

for i = 1:length(cores)
    % read in 14C depths and calibrated ages for each core
    [depth_14c, id1]=read_radiocarbon(cores{i,1});
    [age_14c_2p5, age_14c, age_14c_97p5, id2]=read_14c_calibrated(cores{i,1});
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
    cores_14c{i}=[depth_14c age_14c Bage_14c];
   

    
    %Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
    cores_sr{i}=diff(depth_14c)./diff(Bage_14c);
    cores_nsr(i)=length(cores_sr{i});

	% mean sed rate for core
    cores_meansr(i)=(max(depth_14c)-min(depth_14c))/(max(Bage_14c)-min(Bage_14c));
    
    % Divide by mean sed rate to convert to sed rate ratio
	sed_rate_ratio=cores_sr{i}/cores_meansr(i);
    ind1=find(~isnan(sed_rate_ratio));
    
    % Duration of sed rate
    dur=diff(Bage_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2)];
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
    
    cores_ratio_age{i}=[sed_rate_ratio(ind1) Bage_14c(ind1);...
        sed_rate_ratio(end) Bage_14c(end)];
end

for j = 1:length(cores2)
    [radio_depthage, Bchron_depthage]=func_combine_files_copy(j);
    i=length(cores)+1;
    cores{i,1}=cores2{j};
    %cores_14c{i}=radio_depthage;
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

    cores_14c{i}=[depth_14c age_14c Bage_14c];

    %Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
    cores_sr{i}=diff(depth_14c)./diff(Bage_14c);
    cores_nsr(i)=length(cores_sr{i});
    
	% mean sed rate for core
    cores_meansr(i)=(max(depth_14c)-min(depth_14c))/(max(Bage_14c)-min(Bage_14c));
    
    % Divide by mean sed rate to convert to sed rate ratio
	sed_rate_ratio=cores_sr{i}/cores_meansr(i);
    ind1=find(~isnan(sed_rate_ratio));
    
    % Duration of sed rate
    dur=diff(Bage_14c);
    ind2=find(dur~=0);
    cores_dur{i}=[sed_rate_ratio(ind1) dur(ind2)];
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];
    
    cores_ratio_age{i}=[sed_rate_ratio(ind1) Bage_14c(ind1);...
        sed_rate_ratio(end) Bage_14c(end)];

end

MIN_SPACING=.5; MAX_SPACING=4;
ind=find(cores_meansr>8 & cores_nsr>1);
autocc1=[]; autocc1a=[]; autocc2=[]; autocc2a=[]; 
autocc3=[]; autocc3a=[];
all_sr1=[]; all_sr2=[]; all_dur=[];

for i=1:length(ind)
    depth_age=cores_14c{ind(i)};
    depth=depth_age(:,1);
    age=depth_age(:,3);
    ratio_age=cores_ratio_age{ind(i)};    
    
    asr=cores_dur{ind(i)};
    jnd=find(asr(:,2)>MIN_SPACING & asr(:,2)<MAX_SPACING);
    ratio=ratio_age(:,1);
    ratio=ratio(jnd);
    all_sr1=[all_sr1; ratio];
    all_dur=[all_dur; asr(jnd,2)];
    
    % autocorrelation of ratio between each pair of 14C estimates
    if length(ratio)>3
        cc=corrcoef(ratio(1:end-1),ratio(2:end));
        autocc1(end+1,:)=cc(1,2);
        cc=corrcoef(log(ratio(1:end-1)),log(ratio(2:end)));
        autocc1a(end+1,:)=cc(1,2);
    end
    
    if length(ratio)>4
        cc=corrcoef(ratio(1:end-2),ratio(3:end));
        autocc2(end+1,:)=cc(1,2);
        cc=corrcoef(log(ratio(1:end-2)),log(ratio(3:end)));
        autocc2a(end+1,:)=cc(1,2);
    end
    
    age1=ratio_age(:,2);
    t1=ceil(age1(1));
    t2=floor(age1(end));
    t=[t1:t2];
    
    %eliminate multiple 14c estimates from same depth
    % estimate depth every 1 kyr through core
    [C_uniq,I_data,I_uniq]=unique(age); 
    depth_ev=interp1(age(I_data),depth(I_data),t);
    sr_est=diff(depth_ev);
    mean_sr=(depth_ev(end)-depth_ev(1))/(t2-t1);
    sr_est=sr_est/mean_sr;
    
    % eliminate duplicate ratio estimates due to widely spaced 14C data
    sr_uniq=unique(sr_est);
    all_sr2=[all_sr2; sr_uniq'];
    if length(sr_uniq)>3
        cc=corrcoef(sr_uniq(1:end-1),sr_uniq(2:end));
        autocc3(end+1,:)=[cc(1,2) length(sr_uniq)];
        cc=corrcoef(log(sr_uniq(1:end-1)),log(sr_uniq(2:end)));
        autocc3a(end+1,:)=[cc(1,2) length(sr_uniq)];
    end
end

[mean(autocc1a(:,1)) mean(autocc1(:,1));
    mean(autocc2a(:,1)) mean(autocc2(:,1));
    mean(autocc3a(:,1)) mean(autocc3(:,1))]

[length(all_sr1) length(all_sr2)]
[exp(mean(log(all_sr1))) mean(log(all_sr1)) std(log(all_sr1));
exp(mean(log(all_sr2))) mean(log(all_sr2)) std(log(all_sr2))]

figure(3)
subplot(131)
hist(log(all_sr1),17)
[N,X]=hist(log(all_sr1),17);
hold on
pdf_est=normpdf([-2:.1:2],mean(log(all_sr1)),std(log(all_sr1)));
plot([-2:.1:2],314*pdf_est*mean(diff(X)),'r')
hold off

subplot(132)
hist(log(all_sr2),17)
[N,X]=hist(log(all_sr2),17);
hold on
pdf_est=normpdf([-2:.1:2],mean(log(all_sr2)),std(log(all_sr2)));
plot([-2:.1:2],711*pdf_est*mean(diff(X)),'r')
hold off
%hist(all_dur,17)
subplot(133)
qqplot(log(all_sr2))
%qqplot(log(all_sr1))
