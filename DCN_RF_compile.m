clear

vidFile = '2022_07_05_0020-21_cat.avi';
abfFile = vidFile(1:15);
vidReader = VideoReader(vidFile);
vidPlayer = vision.DeployableVideoPlayer;
results = struct('Boxes',[],'Score',[]);
ctr = [];

data.all = abf2load(strcat(abfFile,'.abf'));
n.epochs = size(data.all,3);
n.frames = vidReader.Duration*vidReader.FrameRate;
% times.FrameStart = n.frames./n.epochs;
times.EpochFrames = 1; % frames per trial
times.EpochStart = (1:times.EpochFrames:n.frames);

times.EpochStart = times.EpochStart(1:n.epochs);

% analyze recorded data

Fs = 20000;

data.all = data.all(:,:,:);

data.EC = shiftdim(data.all(:,1,:),2)';
data.camTrig = shiftdim(data.all(:,2,:),2)';
data.LED = shiftdim(data.all(:,4,:),2)';
data.sol = shiftdim(data.all(:,3,:),2)';
data.cont = shiftdim(data.all(:,4,:),2)';
n.sweeps = size(data.EC,2);

% plot

plot(data.EC(:,18))

%% additional filtering, for units
% 
% d=fdesign.highpass('N,Fc',150,200,20000); % N,6dB,sampling rate
% designmethods(d);
% Hd = design(d);
% 
% data.EC = filter(Hd,data.EC);

% collect spiking
clear sel

thresh = 0.4;

data.dEC = diff(data.EC);

times.stimON = find((diff((data.LED(:,2)>2))>0)==1);
% times.stimOFF = 1.539;
% times.stimON = 1.4542;
% times.stimOFF = 1.47935;
% 
if thresh<0
    spikes = data.EC<thresh;
else
    spikes = data.EC>thresh; 
%     spikes = abs(data.EC)>thresh; 
end

spikes = diff(spikes)>0;
spikes = [zeros(1,size(spikes,2));spikes];

% for ii = 1:n.sweeps
%     times.spikes{ii} = find(diff(spikes(:,ii))>0 == 1);
%     n.spikes.base(ii) = sum(spikes(1:Fs*times.stimON,ii));
%     n.spikes.stim(ii) = sum(spikes(times.stimON*Fs:(times.stimOFF*Fs+Fs*0.05),ii));
%     n.freq.base(ii) = n.spikes.base(ii)./times.stimON;
%     n.freq.stim(ii) = n.spikes.stim(ii)./(times.stimOFF-times.stimON+0.05);
%     n.freq.dF(ii) = n.freq.stim(ii)-n.freq.base(ii);
%     n.spikes.dS(ii) = n.spikes.stim(ii)-n.spikes.base(ii);
% end

trials.LED = data.LED>0.05;
trials.LED = sum(trials.LED);
sel.LED = trials.LED>500;
times.data = (1:length(data.EC))/Fs;

mult = times.stimON;
mult = repmat(mult,1,Fs*0.08);
mult = mult+repmat((1:Fs*0.08),size(mult,1),1);
mult = mult'-0.02*Fs;

temp = spikes(mult,:);

hists.stim = temp;
hists.stim = reshape(temp,Fs*0.08,length(times.stimON),size(data.EC,2));

data.stim = data.EC(mult,:);
data.stim = reshape(data.stim,Fs*0.08,length(times.stimON),size(data.EC,2));
times.datastim = (1:length(data.stim))./Fs-0.02;
% summary.spikes.base = sum(spikes(Fs*(times.stimON-0.1):Fs*(times.stimON),:),1);
temp = reshape(temp,Fs*0.08,length(times.stimON),size(data.EC,2));
summary.spikes.stim = sum(temp((Fs*0.04):(Fs*0.06),:,:),1);
summary.spikes.stim = sum(summary.spikes.stim,2);
summary.spikes.stim = shiftdim(summary.spikes.stim,1);

% mult = times.stimON;
% mult = repmat(mult,1,Fs*0.04);
% mult = mult+repmat((1:Fs*0.04),size(mult,1),1);
% mult = mult';
% mult = mult-0.04*Fs;
% temp = spikes(mult,:);
% 
summary.spikes.base = sum(temp(1:Fs*0.02,:,:),1);
summary.spikes.base = sum(summary.spikes.base,2);
summary.spikes.base = shiftdim(summary.spikes.base,1);

summary.spikes.dF = summary.spikes.stim-summary.spikes.base;

times.discrim = (1:Fs*0.002)-Fs*0.001;
figure, hold on
clear mult
for ii = 1:size(spikes,2)
    times.spikes{ii} = find(spikes(:,ii) == 1);
    mult = times.spikes{ii}';
    mult = repmat(mult,Fs*0.002,1);
    mult = mult+repmat((1:Fs*0.002)',1,size(mult,2));
    mult = mult - Fs*0.001;
    if sum(mult(:))>=0
        mult(mult<=0) = 1;
    end
    discrim.spikes{ii} = data.EC(mult,ii);
    discrim.spikes{ii} = reshape(discrim.spikes{ii},Fs*0.002,length(times.spikes{ii}));
    plot(discrim.spikes{ii});
end     



%% run probe detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,subplot(1,2,1);
set(gcf,'Position',[200 200 1500 600])

% clear results scores bboxes I ctr

for ii = 1:length(times.EpochStart);
    
   
    I = read(vidReader,round(times.EpochStart(ii)));
%   I2 = read(vidReader,round(124));
    subplot(1,2,1),imshow(I)
    subplot(1,2,2),imshow(I2)
    ctr(ii,:) = ginput(1);
    
end
% 
v = read(vidReader,round(times.EpochStart(end-5)));
subplot(1,2,1),hold on,imshow(v),scatter(ctr(:,1),ctr(:,2),9,summary.spikes.stim','filled')
colormap 'spring'




%% trial plotter
figure, set(gcf,'Position',[200 200 1500 600])
times.data = (1:length(data.EC))/Fs;

BWR = [0,0,1;0.0322580635547638,0.0322580635547638,1;0.0645161271095276,0.0645161271095276,1;0.0967741906642914,0.0967741906642914,1;0.129032254219055,0.129032254219055,1;0.161290317773819,0.161290317773819,1;0.193548381328583,0.193548381328583,1;0.225806444883347,0.225806444883347,1;0.258064508438110,0.258064508438110,1;0.290322571992874,0.290322571992874,1;0.322580635547638,0.322580635547638,1;0.354838699102402,0.354838699102402,1;0.387096762657166,0.387096762657166,1;0.419354826211929,0.419354826211929,1;0.451612889766693,0.451612889766693,1;0.483870953321457,0.483870953321457,1;0.516129016876221,0.516129016876221,1;0.548387110233307,0.548387110233307,1;0.580645143985748,0.580645143985748,1;0.612903237342835,0.612903237342835,1;0.645161271095276,0.645161271095276,1;0.677419364452362,0.677419364452362,1;0.709677398204804,0.709677398204804,1;0.741935491561890,0.741935491561890,1;0.774193525314331,0.774193525314331,1;0.806451618671417,0.806451618671417,1;0.838709652423859,0.838709652423859,1;0.870967745780945,0.870967745780945,1;0.903225779533386,0.903225779533386,1;0.935483872890472,0.935483872890472,1;0.967741906642914,0.967741906642914,1;1,1,1;1,0.968750000000000,0.968750000000000;1,0.937500000000000,0.937500000000000;1,0.906250000000000,0.906250000000000;1,0.875000000000000,0.875000000000000;1,0.843750000000000,0.843750000000000;1,0.812500000000000,0.812500000000000;1,0.781250000000000,0.781250000000000;1,0.750000000000000,0.750000000000000;1,0.718750000000000,0.718750000000000;1,0.687500000000000,0.687500000000000;1,0.656250000000000,0.656250000000000;1,0.625000000000000,0.625000000000000;1,0.593750000000000,0.593750000000000;1,0.562500000000000,0.562500000000000;1,0.531250000000000,0.531250000000000;1,0.500000000000000,0.500000000000000;1,0.468750000000000,0.468750000000000;1,0.437500000000000,0.437500000000000;1,0.406250000000000,0.406250000000000;1,0.375000000000000,0.375000000000000;1,0.343750000000000,0.343750000000000;1,0.312500000000000,0.312500000000000;1,0.281250000000000,0.281250000000000;1,0.250000000000000,0.250000000000000;1,0.218750000000000,0.218750000000000;1,0.187500000000000,0.187500000000000;1,0.156250000000000,0.156250000000000;1,0.125000000000000,0.125000000000000;1,0.0937500000000000,0.0937500000000000;1,0.0625000000000000,0.0625000000000000;1,0.0312500000000000,0.0312500000000000;1,0,0];



for f = 1:inf

v = read(vidReader,times.EpochStart(2));
% I = read(vidReader,1);
v = I2;
if ~exist('final')
    subplot(2,2,[1,3]),hold on,imshow(v),scatter(ctr(:,1),ctr(:,2),9,summary.spikes.stim','filled'), %zoom(2)
else
    subplot(2,2,[1,3]),hold on,imshow(v),scatter(final.probe(:,1),final.probe(:,2),9,summary.spikes.stim','filled')
% subplot(1,2,2),hold on,imshow(v),plot(ctr(:,1),ctr(:,2),'.g')
end
set(gca,'colormap',BWR)
caxis([-5 5]);

sel = ginput(1);

if ~exist('final')
d = sqrt((sel(1)-ctr(:,1)).^2+(sel(2)-ctr(:,2)).^2);
else
    d = sqrt((sel(1)-final.probe(:,1)).^2+(sel(2)-final.probe(:,2)).^2);
end

% [A idx] = min(d);
% temp = hists.stim(:,:,idx);
% temp(temp == 0) = nan;
% temp = temp+repmat((1:size(hists.stim,2)),length(temp),1);
% times.hists = (1:length(temp))/Fs
% 
% subplot(2,2,2),plot(times.hists,temp,'.black')
% title(strcat('trial ',num2str(idx)))

[A idx] = min(d);

subplot(2,2,2),cla,plot(times.datastim,data.stim(:,:,idx))
% plot(times.data,data.EC(:,idx),'black')
title(strcat('trial ',num2str(idx),' spk count:', num2str(summary.spikes.stim(idx))))

axis([-0.02 0.04 min(min(data.stim(:,:,idx)))*1.5 max(max(data.stim(:,:,idx)))])

% axis([0 length(temp)/Fs 1 size(temp,2)+2])

end

%% run manual paw detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,subplot(1,2,1);
set(gcf,'Position',[200 200 1500 600])

clear results scores bboxes I 

for ii = 1:n.epochs;
    
    I = read(vidReader,times.EpochStart(ii));
    imshow(I)
    
    paw.ctr(ii,:) = ginput(1);
    % get data

end
%%
final.dpaw = paw.ctr - repmat(nanmean(paw.ctr(6:10,:)),length(paw.ctr),1);
final.dpaw(isnan(final.dpaw)) = 0;
final.probe = ctr - final.dpaw;



