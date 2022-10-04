clear

% concatenate

abffiles.all = getfilenamese('Z:\GintyLab\Turecek\data\2022_08_08_Adv_AcrAcr_ICPN','*.abf');
temp = [];

sweepSel = [];

for f = [1]+1;

    imp = abf2load(abffiles.all{f});
    temp = [temp;shiftdim(imp,2)];
    
end
%%
ff = f;
data.all = shiftdim(temp,1);

if ~isempty(sweepSel)
    sweepSel = sweepSel(sweepSel<size(data.all,3));
    data.all = data.all(:,:,sweepSel);
end
    
Fs = 20000;

data.all = data.all(:,:,:);

data.EC = shiftdim(data.all(:,1,:),2)';
data.camTrig = shiftdim(data.all(:,4,:),2)';
data.LED = shiftdim(data.all(:,4,:),2)';
data.sol = shiftdim(data.all(:,2,:),2)';
n.sweeps = size(data.EC,2);

trials.LED = data.LED>0.02;
trials.LED = sum(trials.LED);
sel.LED = trials.LED>20;
sel.base = trials.LED<20;

data.all = data.all(:,:,or(sel.LED,sel.base));

n.sweeps = size(data.EC,2);

clear temp1 temp2;

%% plot

plot(data.EC(:,27))

%% filtering for units

d=fdesign.highpass('N,Fc',150,200,20000); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);

data.fEC = filter(Hd,data.EC);



%% collect spiking

thresh = -0.3;

data.dEC = diff(data.EC);

temp = 0.9319*Fs;
temp(2) = temp+Fs*0.1;
times.stimON = temp(1)./Fs;
times.stimOFF = temp(end)./Fs; clear temp

if thresh<0
    spikes = data.EC<thresh;
else
    spikes = data.EC>thresh;
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
end


times.data = (1:length(data.EC))/Fs;

%% generate histograms, rasters

wind = Fs*0.001;

clear histog

for f = 1:size(spikes,2)
    temp = reshape(spikes(:,f),[wind,length(spikes)/wind]);
    histog(:,:,f) = temp;
end

histog = sum(histog,1);
histog = shiftdim(histog,1);

times.histog= (1:size(histog,1))/wind;
plot(times.histog,sum(histog,2))

if isfield(sel,'base')
trials.base = ~~sel.base;
else
    trials.base = ~sel.LED;
end
trials.post = ~~sel.LED;

hists.base = sum(histog(:,trials.base),2)./sum(trials.base);%./wind*Fs;
hists.post = sum(histog(:,trials.post),2)./sum(trials.post);%./wind*Fs;


times.hist = ((wind/Fs):(wind/Fs):(length(hists.base)*(wind/Fs)));

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
    
    if trials.base(ii) == 1
        raster.xBase = [raster.xBase; temp];
        raster.yBase = [raster.yBase;ones(length(temp),1)*j];
        j = j+1;
    else
        raster.xPost = [raster.xPost; temp];
        raster.yPost = [raster.yPost;ones(length(temp),1)*k];
        k = k+1;
    end
end

raster.xSort = [raster.xBase;raster.xPost];
raster.ySort = [raster.yBase;raster.yPost+j];



times.disp.raststim.x = double(trials.post);
times.disp.raststim.x(times.disp.raststim.x == 0) = nan;
times.disp.raststim.x = times.disp.raststim.x.*times.stimON-0.02;
times.disp.raststim.y = (1:length(trials.post));
%
figure
set(gcf,'Position',[200 200 1500 600])

ax2 = subplot(5,1,3);
plot(times.hist,hists.base);
hold on
plot(times.hist,hists.post,'r')

ax3 = subplot(5,1,4);
plot(times.data,data.LED)

ax4 = subplot(5,1,5);
subplot(5,1,5),plot(times.data,data.sol(:,1),'color',[0.9 0.9 0.9]);
hold on
times.soladj = Fs*0.00685;
plot(times.data,[zeros(times.soladj,1);data.sol(1:(end-times.soladj),1)],'color',[0.2 0.2 0.2]);

ax1 = subplot(5,1,2);
subplot(5,1,2),plot(raster.xSort,raster.ySort,'.black')
hold on
plot(times.stimON-0.03,j,'>r')
set(gca,'ydir','reverse');
axis([1 1.8 0 ii])

ax0 = subplot(5,1,1);

subplot(5,1,1),plot(raster.x,raster.y,'.black')
set(gca,'ydir','reverse');

axis([0.6 1.7 0 ii])

linkaxes([ax0 ax1 ax2 ax3 ax4],'x');

axis([0.6 1.7 0 ii])

temp = abffiles.all{ff};
temp = strrep(temp,'_',' ');
title(temp)

