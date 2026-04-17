%% Plot figures for each core?
plotCoreFigs = 0; %1 = yes, 0 = no

%% Set Up Cores
%Name cores to be considered
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

cores2={'MD95-2042';'MD01-2416';'MD02-2489';'MD98-2181';'W8709A-13';...
    'MD07-3076';'RC11-83';'GIK17940-2';'V35-5'};

%% Calculate nSR in each core
%Initialize vector to store nSR
all_sed_ratios=[];
cores_14c = cell(1,length(cores));
cores_sr= cell(1,length(cores));
cores_ages= cell(1,length(cores));
cores_ratio_age= cell(1,length(cores));
cores_Bage = cell(1,length(cores));
cores_depth = cell(1,length(cores));

%Cycle through each core
for i = 1:length(cores)
    % read in 14C depths and mean calibrated ages for each core
    [depth_14c, id1]=read_radiocarbon(cores{i,1}); %radiocarbon depths
    [age_14c_2p5, age_14c_med, age_14c_97p5, id2]=read_14c_calibrated(cores{i,1}); % median (50%) calibrated age
    cores_14c{i}=[depth_14c age_14c_med]; %store in cell

    % read in depths and ages from Bchron run (sampled at 1-kyr increments
    % I guess).
    [Bdepth,Bage]=read_14c_ranges(cores{i,1}); %Bdepth is not 14C depths but depths sampled by Bchron, Bage is 50% age at each Bdepth
    % subplot(4,2,2)
    % plot(Bdepth, Bage, 'kx')
    % xlabel("Depth");
    % ylabel("Age");
    % title("step2")

    % Bchron age models frequently don't include first and last 14c dates
    % due to even 1-kyr increments of output
    % Add these ages back in
    if depth_14c(1) < min(Bdepth)
        Bdepth=[depth_14c(1); Bdepth]; %#ok<*AGROW>
        Bage=[age_14c_med(1); Bage];
    end
    if depth_14c(end) > max(Bdepth)
        Bdepth=[Bdepth; depth_14c(end)];
        Bage=[Bage; age_14c_med(end)];
    end

    % subplot(4,2,3)
    % plot(Bdepth, Bage, 'kx')
    %     xlabel("Depth");
    % ylabel("Age");
    % title("step3")
    % 
    %Bchron model ages for depth of each 14c measurement
    % This is taking the median age estimates at approx each
    % 1kyr increment, and then using linear interpolation to get an
    % estimate of the age at each radiocarbon depth. This means that if
    % an age is an outlier (and ignored by Bchron) then an
    % age at that depth (that has nothing to do with the radiocarbon date
    % at that depth) will still be read and treated as if it is from that depth.
    Bage_14c=interp1(Bdepth,Bage,depth_14c);
    cores_Bage{i}=Bage_14c;
    cores_depth{i}=depth_14c;

    % Neighbouring Age Pairs
    NP_age = [Bage_14c(1:end-1), Bage_14c(2:end)];
    NP_depth = [depth_14c(1:end-1), depth_14c(2:end)];
    NP_agediff = diff(NP_age, [], 2);
    NPf_ind = NP_agediff >= 0.5 & NP_agediff <= 4;

    NPf_age = NP_age(NPf_ind,:);
    NPf_depth = NP_depth(NPf_ind,:);
    
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
    core_SRandDur = [sed_rate_ratio(ind1) dur(ind2) depint(ind2)];
    cores_dur{i}= core_SRandDur;
    mindur(i)=min(dur(ind2));
    all_sed_ratios=[all_sed_ratios; cores_dur{i}];

    %Plotting histogram
    MIN_SPACING=.5; MAX_SPACING=4;

    indFilt=find(core_SRandDur(:,2)>=MIN_SPACING & core_SRandDur(:,2)<=MAX_SPACING);
    core_SRandDurFilt=core_SRandDur(indFilt,:);

    coreSRtW = [];
    coreSRdW = [];
    tduration=round(10*core_SRandDurFilt(:,2));
    depthint = round(core_SRandDurFilt(:,3));
    for k=1:length(core_SRandDurFilt(:,1))
        if tduration(k) > 0
            coreSRtW = [coreSRtW; core_SRandDurFilt(k,1).*ones(tduration(k), 1)];
        end
        if tduration(k) >0
            coreSRdW = [coreSRdW; core_SRandDurFilt(k,1).*ones(depthint(k), 1)];
        end
    end


    if plotCoreFigs
    figure;
    subplot(4,2,[1 3])
    hold on
    plot(depth_14c, age_14c_med, 'r.', "HandleVisibility", "off");
    errorbar(depth_14c, age_14c_med, age_14c_med-age_14c_2p5, age_14c_97p5-age_14c_med, 'vertical', 'Color', 'r', 'LineStyle','none', "DisplayName", "Calibrated Age")
    xlabel("Depth");    
    ylabel("Age");
    ylims = ylim;
    box on;

    hold on
    plot(depth_14c, Bage_14c, 'k.', "DisplayName", "Age Used")
    xlabel("Depth");
    ylabel("Age");
    legend("Location", "best")

    subplot(4,2,[2 4])
    plot(depth_14c, Bage_14c, 'k.', "DisplayName","Age Used")
    xlabel("Depth");
    ylabel("Age");
    ylim(ylims);
    box on;
    hold on
    plot(NPf_depth', NPf_age', '-', 'LineWidth', 1,"HandleVisibility", "off")
    plot([min(depth_14c) max(depth_14c)], [min(Bage_14c) max(Bage_14c)], 'k--', 'LineWidth', 1,"DisplayName", "Overall SR")
    legend("Location", "best")

    subplot(4,2,[5 6]);
    histogram(coreSRtW, 0:0.1:10);
    xlim([0 6])
    xlabel("nSR")
    ylabel("0.1kyr increments")

    subplot(4,2,[7 8]);
    histogram(coreSRdW, 0:0.1:10);
    xlim([0 6])
    xlabel("nSR")
    ylabel("1cm increments")
    
    sgtitle(string(cores{i}))
    end
end

%Do the same thing for cores in cores2
for j = 1:length(cores2)
    %Get depth and median calibrated ages of radiocarbon depths, and depth
    %and median ages of 1kyr increments in Bchron
    [radio_depthage, Bchron_depthage]=func_combine_files_copy(j); %combine info from multiple files (uncertain why some cores were run in by splitting up and running multiple Bchron age models)
    i=length(cores)+1;
    cores{i,1}=cores2{j};
    cores_14c{i}=radio_depthage;
    depth_14c=radio_depthage(:,1);
    age_14c_med=radio_depthage(:,2);
    age_14c_2p5 = radio_depthage(:,3);
    age_14c_97p5 = radio_depthage(:,4);
    Bdepth=Bchron_depthage(:,1);
    Bage=Bchron_depthage(:,2);

    
    
    % Bchron age models frequently don't include first and last 14c dates
    % due to even 1-kyr increments of output
    % Add these ages back in
    if depth_14c(1) < min(Bdepth)
        Bdepth=[depth_14c(1); Bdepth];
        Bage=[age_14c_med(1); Bage];
    end
    if depth_14c(end) > max(Bdepth)
        Bdepth=[Bdepth; depth_14c(end)];
        Bage=[Bage; age_14c_med(end)];
    end
    
    %Bchron model ages for depth of each 14c measurement
    % This is taking the median age estimates at approx each
    % 1kyr increment, and then using linear interpolation to get an
    % estimate of the age at each radiocarbon depth...?! This means that if
    % an age is an outlier (and ignored by Bchron) then this will get an
    % age at that depth (that has nothing to do with the radiocarbon date
    % at that depth) and treat it as if it is from that depth.
    Bage_14c=interp1(Bdepth,Bage,depth_14c);
    cores_Bage{i}=Bage_14c;
    cores_depth{i}=depth_14c;

    

    % Neighbouring Age Pairs
    NP_age = [Bage_14c(1:end-1), Bage_14c(2:end)];
    NP_depth = [depth_14c(1:end-1), depth_14c(2:end)];
    NP_agediff = diff(NP_age, [], 2);
    NPf_ind = NP_agediff >= 0.5 & NP_agediff <= 4;

    NPf_age = NP_age(NPf_ind,:);
    NPf_depth = NP_depth(NPf_ind,:);

    

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

    if plotCoreFigs
        figure;
    subplot(4,2,[1 3])
    hold on
    plot(depth_14c, age_14c_med, 'r.', "HandleVisibility", "off");
    errorbar(depth_14c, age_14c_med, age_14c_med-age_14c_2p5, age_14c_97p5-age_14c_med, 'vertical', 'Color', 'r', 'LineStyle','none', "DisplayName", "Calibrated Age")
    xlabel("Depth");    
    ylabel("Age");
    ylims = ylim;
    box on;
    hold on
    plot(depth_14c, Bage_14c, 'k.', "DisplayName", "Age Used")
    xlabel("Depth");
    ylabel("Age");
    legend("Location", "best")

    subplot(4,2,[2 4])
    plot(depth_14c, Bage_14c, 'k.', "DisplayName","Age Used")
    xlabel("Depth");
    ylabel("Age");
    ylim(ylims);
    box on;
    hold on
    plot(NPf_depth', NPf_age', '-', 'LineWidth', 1,"HandleVisibility", "off")
    plot([min(depth_14c) max(depth_14c)], [min(Bage_14c) max(Bage_14c)], 'k--', 'LineWidth', 1,"DisplayName", "Overall SR")
    legend("Location", "best")
    
      subplot(4,2,[5 6]);
    histogram(coreSRtW, 0:0.1:10);
    xlim([0 6])
    xlabel("nSR")
    ylabel("0.1kyr increments")


    subplot(4,2,[7 8]);
    histogram(coreSRdW, 0:0.1:10);
    xlim([0 6])
    xlabel("nSR")
    ylabel("1cm increments")

    
    sgtitle(string(cores{i}))
    end
end

%% Find cores with a mean SR >=8 and where minimum spacing <=4kyr
asr=[];
indcores=find(cores_meansr>=8 & mindur<=4);
for i=1:length(indcores)
    asr=[asr; cores_dur{indcores(i)}; [NaN NaN NaN]];
end

asrnan=asr;

MIN_SPACING=.5; MAX_SPACING=4;
ind=find(asr(:,2)>=MIN_SPACING & asr(:,2)<=MAX_SPACING);
asr=asr(ind,:);

ind=find((asrnan(:,2)>=MIN_SPACING & asrnan(:,2)<=MAX_SPACING) | isnan(asrnan(:,1)) );
asrnan=asrnan(ind,:);


disp('----------------------------------------------------')
disp(['MIN SPACING = ' num2str(MIN_SPACING)])
disp(['MAX SPACING = ' num2str(MAX_SPACING)])
disp(['Number of cores = ' num2str(length(indcores))])

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
mean_srw=mean(srw);
mean_srdw=mean(srdw);

%Save dts
save("Lin2014_dts.mat", ["asr"], '-mat');

figure(2)
subplot(211)    
[n, ~, ~]=histcounts(srw,ratio_bins);
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
subplot(2,1,1)
binWidth = 0.065;
middlebin = [0-binWidth/2 0+binWidth/2];
upperbins = middlebin(2)+binWidth*(1:35);
lowerbins = middlebin(1)+binWidth*(-40:-1);
binEdges2 = [lowerbins, middlebin, upperbins];
histogram(log(srw),binEdges2)
set(gca,'XTick',[-1.38 -0.0812 0 0.0816 1.38],'XTickLabel',[.25 .92 1 1.08 4])
ax = gca;                % get current axes
ax.XAxis.FontSize = 6;  
ax.YAxis.FontSize = 6;  
axis tight
xlabel('log(nSR)')
ylabel('0.1 kyr increments')
%title('Reproducing nSR histogram from Lin et al')

%%
subplot(2,1,2)
histogram(log(srdw),[-2.242:.065:2])
%set(gca,'XTick',[-1.38 -0.0812 0 0.0816 1.38],'XTickLabel',[.25 .92 1 1.08 4])
axis tight
xlabel('log(nSR)')
ylabel('1cm counts')
title('Reproducing nSR histogram from Lee et al')

