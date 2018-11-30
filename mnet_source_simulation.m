function mnet_source_simulation(cfg, data, sourcemodelfile, headmodelfile, transformfile)
%% Simulate Signal
load(sourcemodelfile);
load(headmodelfile);
use_my_library('spm',0);
use_my_library('ft',1);
cfg.timepoint = 0:1/cfg.fsample:cfg.time;

signal = zeros(size(cfg.timepoint));
signal(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition))) = ... 
    (0.5 - 0.5*cos(2*pi*(cfg.timepoint(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition)))-cfg.wave_start)/cfg.wavelet_length)) ... 
    .*(cfg.wave_amp*sin(2*pi*cfg.wave_freq*(cfg.timepoint(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition)))-cfg.wave_start)));
signal = signal +2*(rand(size(cfg.timepoint))-0.5)*cfg.noise_amp;
figure;
subplot(2,3,1);
plot(cfg.timepoint,signal);
title('Simulted Signal');
%% Show Simulated Position on Cortical Mesh
use_my_library('ft',0);
use_my_library('spm',1);

fid = fopen(transformfile);
eval(fread(fid,[1 inf], '*char'));
fclose(fid);
clear fid

fromMNI = transform.spm2bti;
roi = cfg.ROI;
roi(4,:) = 1;
roi = fromMNI*roi;
roi = roi(1:3,:);

for i =1:size(roi,2)
    tempDistance = [];
    tempDistance = sourcemodel2d.pos*10 - repmat(roi(:,i)',size(sourcemodel2d.pos,1),1);
    tempDistance = sum(tempDistance.^2,2);
    roiInd(i) = find(min(tempDistance)==tempDistance);
    neighbors{i} = find_neighbors_jhs(roiInd(i), sourcemodel2d.tri,2);
end
clear tempDistance i

use_my_library('spm',0);
use_my_library('ft',1);

bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;

source_sim = zeros(size(bnd.pnt,1),1);
for i = 1:size(roi,2)
    source_sim(neighbors{i}) = sum(abs(fft(signal)))*i; % multiply i for just separate color for each roi 
end
mask = ~~(source_sim);
alpha = 0.9;
mask = mask + (~mask)*alpha;
subplot(1,3,2);
ft_plot_mesh(bnd);
ft_plot_mesh(bnd, 'vertexcolor', source_sim, 'facealpha',mask,'colormap','jet');
lighting gouraud
camlight

title('Simulated Source Location');
%% Simulate Channel Signals
tmpcfg = [];
tmpcfg.headmodel = headmodel;
tmpcfg.grid = sourcemodel2d;
tmpcfg.grad = data.grad;
tmpcfg.channel = data.label;

leadfield = ft_prepare_leadfield(tmpcfg);

% set simultate data structure.
data_sim = data;
data_sim.fsample = cfg.fsample;
data_sim.time = [];
for i = 1:cfg.ntrial
    data_sim.time{i} = cfg.timepoint;
end

use_my_library('ft',0);
use_my_library('spm',1);
for i=1:length(roiInd)
    neighbors_1 = setdiff(find_neighbors_jhs(roiInd(i), sourcemodel2d.tri,1),roiInd(i));
    if isempty(cfg.ROI_ori) || (numel(cfg.ROI_ori)~=1 && numel(cfg.ROI_ori)~=size(roi,2))
        source_ori{i} =(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(roiInd(i),:))/norm(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(roiInd(i),:));
    else
        %This part is not yet implemented.
    end
    current_mom{i}(1,:) = source_ori{i}(1).*signal;
    current_mom{i}(2,:) = source_ori{i}(2).*signal;
    current_mom{i}(3,:) = source_ori{i}(3).*signal;
end
for i=1:cfg.ntrial
    chan_signal{i} = zeros(length(data.label),length(cfg.timepoint));
    for j =1:length(roiInd)
        for k = 1:length(neighbors{j})
            chan_signal{i} = chan_signal{i} + leadfield.leadfield{neighbors{j}(k)}*current_mom{j}; %% just sum of signal. Should I use L2-norm?
        end
    end
end
data_sim.trial = chan_signal; %% Simulated channel data.

use_my_library('spm',0);
use_my_library('ft',1);

subplot(2,3,4);
plot(data_sim.time{1}, data_sim.trial{1});

clear i j current_mom neighbors neighbors_1 source_mom source_ori source_sim
clear amplitude chan_signal k roi signal tmpcfg

%% Invert and Display
% Invert with LCMV
tmpcfg                   = [];
tmpcfg.method            = 'lcmv';
tmpcfg.grid              = leadfield;
tmpcfg.headmodel         = headmodel;
tmpcfg.lcmv.lambda = '10%';
tmpcfg.lcmv.fixedori = 'yes';
tmpcfg.lcmv.projectnoise = 'yes';
tmpcfg.lcmv.keepleadfield = 'yes';
source_lcmv = ft_sourceanalysis(tmpcfg, data_sim);
source_lcmv = ft_sourcedescriptives([], source_lcmv);

scale_factor = median(floor(log10(source_lcmv.avg.nai(find(~isnan(source_lcmv.avg.nai))))));
source_lcmv.avg.nai = 10.^(log10(source_lcmv.avg.nai)-scale_factor);

% plot the neural activity index (power/noise)
subplot(2,3,3);
ft_plot_mesh(bnd);
ft_plot_mesh(bnd, 'vertexcolor', source_lcmv.avg.nai, 'facealpha',source_lcmv.avg.nai,'alphalim',[3 8],'colormap','jet');
lighting gouraud
camlight
title(sprintf('Source Reconstruction of Simulated Signal\n(Beamformer LCMV)'));

% Invert with loreta
tmpcfg                   = [];
tmpcfg.method            = 'eloreta';
tmpcfg.grid              = leadfield;
tmpcfg.headmodel         = headmodel;
tmpcfg.eloreta.lambda = 10;
tmpcfg.eloreta.fixedori = 'yes';
tmpcfg.eloreta.projectnoise = 'yes';
tmpcfg.eloreta.keepleadfield = 'yes';
source_eloreta = ft_sourceanalysis(tmpcfg, data_sim);
source_eloreta = ft_sourcedescriptives([], source_eloreta);

% plot the neural activity index (power/noise)
subplot(2,3,6);
ft_plot_mesh(bnd);
ft_plot_mesh(bnd, 'vertexcolor', source_eloreta.avg.pow, 'facealpha',source_eloreta.avg.pow,'alphalim',[0 max(source_eloreta.avg.pow)],'colormap','jet');
lighting gouraud
camlight
title(sprintf('Source Reconstruction of Simulated Signal\n(eLORETA)'));

