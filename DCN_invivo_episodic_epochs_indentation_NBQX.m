clear

% concatenate

abffiles.all = getfilenamese('Z:\GintyLab\Turecek\data\2022_01_24_DCN_NBQX','*.abf');
temp = [];

sweepSel = [];
k = 0;

for f = [20:30]+1;

    imp = abf2load(abffiles.all{f});
    data.sweeps = 18;
    imp = imp(:,:,1:18);
    temp = [temp;shiftdim(imp,2)];
    k = k+1;
    
end
ff = f;

data.all = shiftdim(temp,1);
data.files = k;

if ~isempty(sweepSel)
    sweepSel = sweepSel(sweepSel<size(data.all,3));
    data.all = data.all(:,:,sweepSel);
end

Fs = 20000;

data.all = data.all(:,:,:);

data.EC = shiftdim(data.all(:,1,:),2)';
data.camTrig = shiftdim(data.all(:,2,:),2)';
data.LED = shiftdim(data.all(:,3,:),2)';
data.sol = shiftdim(data.all(:,4,:),2)';

if max(data.sol(:,1)<3.21)
    data.sol = data.sol+0.0513;
end

data.sorted = [];
for f = 0:(data.sweeps-1)
    mult = (1:(data.sweeps):(data.files)*data.sweeps)+f;
    data.sorted = [data.sorted,data.EC(:,mult)];
end

n.sweeps = size(data.EC,2);

% clear artifacts from indenter
data.sorted(16469:16477,:) = NaN;
data.sorted(10469:10477,:) = NaN;


%% plot

plot(data.EC(:,18))

%% filtering for units

d=fdesign.highpass('N,Fc',150,200,20000); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);
backup = data.EC;
data.EC = filter(Hd,data.EC);
data.EC = [data.EC(75:end,:);zeros(75,size(data.EC,2))];

%% filtering for fields

d=fdesign.lowpass('N,Fc',1000,200,20000); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);

raw = filter(Hd,data.EC);



%% collect spiking

thresh = -0.4;

data.dsorted = diff(data.sorted);


temp = find((data.sol(:,1)>0.1) == 1);

times.stimON = temp(1)./Fs;
times.stimOFF = temp(end)./Fs; clear temp

if thresh<0
    spikes = data.sorted<thresh;
else
    spikes = data.sorted>thresh;
end

spikes = diff(spikes)>0;
spikes = [zeros(1,size(spikes,2));spikes];

for ii = 1:n.sweeps
    times.spikes{ii} = find(diff(spikes(:,ii))>0 == 1);
    n.spikes.base(ii) = sum(spikes(1:Fs*times.stimON,ii));
    n.spikes.stim(ii) = sum(spikes(times.stimON*Fs:(times.stimOFF*Fs+Fs*0.05),ii));
    n.freq.base(ii) = n.spikes.base(ii)./times.stimON;
    n.freq.stim(ii) = n.spikes.stim(ii)./(times.stimOFF-times.stimON+0.05);
    n.freq.dF(ii) = n.freq.stim(ii)-n.freq.base(ii);
    n.spikes.dS(ii) = n.spikes.stim(ii)-n.spikes.base(ii);
    
    data.stimAmp(ii) = nanmean(data.sol((times.stimON*Fs+Fs*0.01):(times.stimON*Fs+Fs*0.02),ii),1);
end

data.stimAmp = data.stimAmp(1:data.sweeps);



times.discrim = (1:Fs*0.002)-Fs*0.001;
figure, hold on
for ii = 1:size(spikes,2)
    times.spikes{ii} = find(spikes(:,ii) == 1);
    mult = times.spikes{ii}';
    mult = repmat(mult,Fs*0.002,1);
    mult = mult+repmat((1:Fs*0.002)',1,size(mult,2));
    mult = mult - Fs*0.001;
    mult(mult<1) = 1;
    mult(mult>length(spikes)) = length(spikes);
    discrim.spikes{ii} = data.sorted(mult,ii);
    discrim.spikes{ii} = reshape(discrim.spikes{ii},Fs*0.002,length(times.spikes{ii}));
    plot(discrim.spikes{ii});
end   

dcm.vcomm = [3.24810882100000;3.29823052700000;3.34803946600000;3.39791672700000;3.44846203600000;3.49829526700000;3.54833346700000;3.59864647700000;3.64832939800000;3.69832053100000;3.74868212600000;3.79848802900000;3.84863858200000;3.89880128200000;3.94896398100000;3.99912668100000;4.04928938100000;4.09945208000000;4.14961478000000;4.19977748000000];

data.force = [0,1.10000000000000,2.75000000000000,5.10000000000000,8.75000000000000,13.9750000000000,21.7500000000000,32.5750000000000,46.7500000000000,65.2000000000000,87.7000000000000,114.200000000000,143.375000000000,175.375000000000,207.850000000000,239.900000000000,272.800000000000,300.575000000000];

trials.LED = data.LED>0.1;
trials.LED = sum(trials.LED);
sel.LED = trials.LED>100;
times.data = (1:length(data.EC))/Fs;


%% generate histograms, rasters

wind = Fs*0.01;

clear histog hists

for f = 1:size(spikes,2)
    temp = reshape(spikes(:,f),[wind,length(spikes)/wind]);
    histog(:,:,f) = temp;
end

histog = sum(histog,1);
histog = shiftdim(histog,1);

times.histog = (1:size(histog,1))/wind;
% plot(times.histog,sum(histog,2))

trials.base = ~sel.LED;
trials.post = ~~sel.LED;

for f = 1:data.sweeps
    temp = sum(spikes(:,(1:data.files)+(data.files*(f-1))),2);
    mult = repmat((1:wind)',1,length(spikes)/wind);
    mult = mult + repmat(0:wind:(wind*(size(mult,2)-1)),size(mult,1),1);
    hists.ramp(:,f) = sum(temp(mult),1)./data.files;
end


times.hist = ((wind/Fs):(wind/Fs):(length(hists.ramp)*(wind/Fs)));

% generate rasters

raster.x = [];
raster.y = [];
raster.xBase = [];
raster.yBase = [];
raster.xPost = [];
raster.yPost = [];
raster.xSort = [];
raster.ySort = [];
j = 1;
k = 1;

for ii = 1:size(data.EC,2);
        
    temp = find(spikes(:,ii) == 1)./Fs;  
    raster.x = [raster.x;temp];
    raster.y = [raster.y;ones(length(temp),1)*ii];

end



raster.xSort = [raster.xBase;raster.xPost];
raster.ySort = [raster.yBase;raster.yPost+j];



times.disp.raststim.x = double(trials.post);
times.disp.raststim.x(times.disp.raststim.x == 0) = nan;
times.disp.raststim.x = times.disp.raststim.x.*times.stimON-0.02;
times.disp.raststim.y = (1:length(trials.post));


mult = (2:2:data.sweeps);
figure
set(gcf,'Position',[500 000 600 800])

ax1 = subplot(10,2,[2 4 6 8 10 12 14 16 18]);
hold on
k = 0;
clear temp
for f = mult
    plot(times.hist-times.stimON, hists.ramp(:,f)-k,'black')
%     plot(times.hist-times.stimON,nanmean([hists.ramp(:,f),hists.ramp(:,f-1)],2)-k,'black')
    k = k+1;
    temp{k} = num2str(round(data.force(mult(k))));
end
temp = fliplr(temp);
yticks((-k+1):0)
yticklabels(temp)
ylim([-k+0.5 1])

ax2 = subplot(10,2,[1 3 5 7 9 11 13 15 17]);
hold on
temp = data.files.*(1:data.sweeps);
plot(repmat([-1 5],data.sweeps,1)',[temp;temp],'--','color',[0.95 0.95 0.95])
plot(raster.x-times.stimON,raster.y,'.black'),set(gca,'ydir','reverse'),ylim([0 max(raster.y)+1])

yticks(temp)
for ii = 1:length(temp)
    lab{ii} = num2str(round(data.force(ii),1));
end
yticklabels(lab)
% ylim([-k+0.5 1])


% plot(
temp = abffiles.all{ff};
temp = strrep(temp,'_',' ');
title(temp)

ax3 = subplot(10,2,19);
plot(times.data-times.stimON,data.sol(:,1),'black')
ax4 = subplot(10,2,20);
plot(times.data-times.stimON,data.sol(:,1),'black')

linkaxes([ax1 ax2 ax3 ax4],'x');
xlim([-0.01 0.05])
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

mult = (2:2:data.sweeps);
figure
set(gcf,'Position',[500 000 600 800])

ax1 = subplot(10,2,[2 4 6 8 10 12 14 16 18]);
hold on
k = 0;
clear temp
for f = mult
%     N = matlab.lang.makeValidName(strcat('ax',num2str(f)));
    %     N = genvarname(strcat('ax',num2str(f)));
    
%     eval([N ' = subplot(max(mult),2,f)'])
    plot(times.hist-times.stimON, hists.ramp(:,f)-k,'black')
%     plot(times.hist-times.stimON,nanmean([hists.ramp(:,f),hists.ramp(:,f-1)],2)-k,'black')
    k = k+1;
    temp{k} = num2str(round(data.force(mult(k))));
end
temp = fliplr(temp);
yticks((-k+1):0)
yticklabels(temp)
ylim([-k+0.5 1])

ax2 = subplot(10,2,[1 3 5 7 9 11 13 15 17]);
hold on
temp = data.files.*(1:data.sweeps);
plot(repmat([-1 5],data.sweeps,1)',[temp;temp],'--','color',[0.95 0.95 0.95])
yticks(temp)
for ii = 1:length(temp)
    lab{ii} = num2str(round(data.force(ii),1));
end
yticklabels(lab)

plot(raster.x-times.stimON,raster.y,'.black'),set(gca,'ydir','reverse'),ylim([0 max(raster.y)+1])
temp = abffiles.all{ff};
temp = strrep(temp,'_',' ');
title(temp)



ax3 = subplot(10,2,19);
plot(times.data-times.stimON,data.sol(:,1),'black')
ax4 = subplot(10,2,20);
plot(times.data-times.stimON,data.sol(:,1),'black')

linkaxes([ax1 ax2 ax3 ax4],'x');
xlim([-0.1 0.5])


clear i j  N ans temp k f lab ii i j k
