clear

% concatenate

abffiles.all = getfilenamese('Z:\GintyLab\Turecek\data\2022_06_09_ICPN_practice','*.abf');
temp = [];

sweepSel = [];
k = 0;

for f = [1]+1;

    imp = abf2load(abffiles.all{f});
    temp = [temp;shiftdim(imp,2)];
    k = k+1;
    data.nsweeps(k) = size(imp,3);
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

if size(imp,2)>2

data.EC = shiftdim(data.all(:,1,:),2)';
data.camTrig = shiftdim(data.all(:,2,:),2)';
data.LED = shiftdim(data.all(:,2,:),2)';
data.sol = shiftdim(data.all(:,3,:),2)';
else
    data.EC = shiftdim(data.all(:,1,:),2)';
    data.sol = shiftdim(data.all(:,2,:),2)';
end

for ii = 1:size(data.EC,2)
temp = data.sol(:,ii)>0.1;
temp = find((diff(temp)>0) == 1);
freq(ii) = nanmean(diff(temp));
end

freq = round(Fs./freq);

n.sweeps = size(data.EC,2);
data.sweeps = n.sweeps;
%% plot

plot(data.EC(:,50))

%% filtering for units

d=fdesign.highpass('N,Fc',150,200,20000); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);

data.EC = filter(Hd,data.EC);

%% filtering for fields

d=fdesign.lowpass('N,Fc',1000,200,20000); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);

raw = filter(Hd,data.EC);



%% collect spiking

thresh = -0.3;

data.dEC = diff(data.EC);


% temp = find((data.sol(:,1)>0.1) == 1);
temp = 14620;

times.stimON = temp(1)./Fs;
times.stimOFF = temp(1)+0.1; clear temp



if thresh<0
    spikes = data.EC<thresh;
else
    spikes = data.EC>thresh;
end

spikes = diff(spikes)>0;
spikes = [zeros(1,size(spikes,2));spikes];

for ii = 1:n.sweeps 
    data.stimAmp(ii) = max(data.sol((times.stimON*Fs):(times.stimON*Fs+Fs*0.05),ii),[],1);
end

data.stimAmp = data.stimAmp(1:n.sweeps);

dcm.force = [0;0.257142857000000;1.02857142900000;2.01428571400000;3.25714285700000;5.02857142900000;7.45714285700000;10.6000000000000;15.1000000000000;20.6571428600000;27.9714285700000;35.3714285700000;44.9428571400000;54.5000000000000;65.2285714300000;75.9428571400000;86.1428571400000;96.2857142900000;104.300000000000;103.857142900000]
dcm.vcomm = [3.24810882100000;3.29823052700000;3.34803946600000;3.39791672700000;3.44846203600000;3.49829526700000;3.54833346700000;3.59864647700000;3.64832939800000;3.69832053100000;3.74868212600000;3.79848802900000;3.84863858200000;3.89880128200000;3.94896398100000;3.99912668100000;4.04928938100000;4.09945208000000;4.14961478000000;4.19977748000000];

for ii = 1:length(data.stimAmp)
    temp = abs(data.stimAmp(ii)-dcm.vcomm);
    [A t] = min(temp);
    data.force(ii) = dcm.force(t);
end

times.data = (1:length(data.EC))/Fs;

%% generate histograms, rasters

wind = Fs*0.001;

clear histog hists

for f = 1:size(spikes,2)
    temp = reshape(spikes(:,f),[wind,length(spikes)/wind]);
    histog(:,:,f) = temp;
end

histog = sum(histog,1);
histog = shiftdim(histog,1);

for f = 1:data.files;
    if f == 1
        hists(:,f) = sum(histog(:,1:sum(data.nsweeps(1:f))),2)./sum(data.nsweeps(1:f));
    else
        temp = (sum(data.nsweeps(1:(f-1)))+1):(sum(data.nsweeps(1:f)));
    hists(:,f) = sum(histog(:,temp),2)./length(temp);
    end
    
end



times.histog = (1:size(histog,1))/wind;
plot(times.histog,sum(histog,2))

times.hist = ((wind/Fs):(wind/Fs):(length(hists)*(wind/Fs)));

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


figure('units','normalized','outerposition',[0 0.3 1 0.5])
clear temp

temp = unique(freq);
times.data = (1:length(data.EC))/Fs;


for ff = 1:data.files
    
    
% subplot(1,length(data.files),ff)
axes('Position',[ff/(data.files+1)-0.03 0.67 1/(data.files+1) 0.32])
plot(raster.x,raster.y,'.black')


if ff == 1
axis([times.stimON-0.05 times.stimON+0.15 1 data.nsweeps(ff)])
else
    axis([times.stimON-0.05 times.stimON+0.15 sum(data.nsweeps(1:(ff-1)))+1 sum(data.nsweeps(1:ff))])
end

axes('Position',[ff/(data.files+1)-0.03 0.34 1/(data.files+1) 0.32])
plot(times.data,data.sol(:,sum(data.nsweeps(1:ff))-5),'black')

axis([times.stimON-0.05 times.stimON+0.15 0 4])


axes('Position',[ff/(data.files+1)-0.03 0.05 1/(data.files+1) 0.32])
plot(times.hist,hists(:,ff),'black')

axis([times.stimON-0.05 times.stimON+0.15 0 1.1])


end








