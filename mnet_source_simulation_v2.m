function mnet_source_simulation(cfg, data, sourcemodel, headmodel, transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG Source Localization Code                                            %
%     cfg.time           - Length of simulation                           %
%         fsample        - Sampling frequency for simulation              %
%         wave_start     - The time that sine wave starts.                %
%         wave_freq      - Sine wave frequency                            %
%         wave_amp       - Sime wave amplitude                            %
%         wavelet_length - Wavelet wavelength                             %
%         wavelet_repetition - Repetition of wavelet                      %
%         noise_amp      - Noise amplitude                                %
%         ROI            - Each columns are mni coordinates of each ROI.  %
%                          [3 x number_of_ROIs] matrix                    %
%         ROIname        - Names for each ROI. Cell array                 %
%                          {1 x number_of_ROIs}                           %
%         ntrial         - Number of trials to generate                   %
%         ROI_ori        - Orientation of diples                          %
%                       Default : current randomly orients neighbor vertex%
%     data        - Use channel informations of data
%     sourcemodel - sourcemodel(2D surface model)
%     headmodel   - Volume conduction model
%     transform   - If sourcemodel does not MNI coordinates, transform 
%                   matrix must be given
%                   [4 x 4] matrix to MNI coordinates to sourcemodel
%                   coordinates such that 
%                   sourceCoord = transform * MNICoord;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.09.06 16:22 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate Signal
cfg.timepoint = 0:1/cfg.fsample:cfg.time;

signal = zeros(size(cfg.timepoint));
signal(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition))) = ... 
    (0.5 - 0.5*cos(2*pi*(cfg.timepoint(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition)))-cfg.wave_start)/cfg.wavelet_length)) ... 
    .*(cfg.wave_amp*sin(2*pi*cfg.wave_freq*(cfg.timepoint(find((cfg.timepoint>cfg.wave_start).* (cfg.timepoint<cfg.wave_start+cfg.wavelet_length*cfg.wavelet_repetition)))-cfg.wave_start)));
signal = signal +2*(rand(size(cfg.timepoint))-0.5)*cfg.noise_amp;
figure;
subplot(2,3,1);
plot(cfg.timepoint,signal);
title('Simululated Signal');
%% Show Simulated Position on Cortical Mesh
roi = cfg.ROI;
roi(4,:) = 1;
roi = transform*roi;
roi = roi(1:3,:);
for i =1:size(roi,2)
    tempDistance = [];
    tempDistance = sourcemodel.pos - repmat(roi(:,i)',size(sourcemodel.pos,1),1);
    tempDistance = sum(tempDistance.^2,2);
    roiInd(i) = find(min(tempDistance)==tempDistance);
    neighbors{i} = find_neighbors_jhs(roiInd(i), sourcemodel.tri,6);
end
clear tempDistance i

bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;

source_sim = zeros(size(bnd.pnt,1),1);
for i = 1:size(roi,2)
    source_sim(neighbors{i}) = sum(abs(fft(signal)))*i; % multiply i for just separate color for each roi 
end
mask = ~~(source_sim);
alpha = 0.9;
mask = mask + (~mask)*alpha;
subplot(2,3,2);
ft_plot_mesh(bnd);
ft_plot_mesh(bnd, 'vertexcolor', source_sim, 'facealpha',mask,'colormap','jet');
lighting gouraud
camlight

title('Simulated Source Location');
%% Simulate Channel Signals
tmpcfg = [];
tmpcfg.headmodel = headmodel;
tmpcfg.grid = sourcemodel;
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

for i=1:length(roiInd)
    neighbors_1 = setdiff(find_neighbors_jhs(roiInd(i), sourcemodel.tri,1),roiInd(i));
    if isempty(cfg.ROI_ori) || (numel(cfg.ROI_ori)~=1 && numel(cfg.ROI_ori)~=size(roi,2))
        source_ori{i} = (sourcemodel.pos(neighbors_1(1),:) - sourcemodel.pos(roiInd(i),:))/norm(sourcemodel.pos(neighbors_1(1),:) - sourcemodel.pos(roiInd(i),:));
    elseif numel(cfg.ROI_ori) == 1
        source_ori{i} = cfg.ROI_ori{1};
    elseif numel(cfg.ROI_ori)==size(roi,2)
        source_ori{i} = cfg.ROI_ori{i};
    end
    current_mom{i}(1,:) = source_ori{i}(1).*signal;
    current_mom{i}(2,:) = source_ori{i}(2).*signal;
    current_mom{i}(3,:) = source_ori{i}(3).*signal;
end
for i=1:cfg.ntrial
    chan_signal{i} = zeros(length(data.label),length(cfg.timepoint));
    for j =1:length(roiInd)
        for k = 1:length(neighbors{j})
            if ~isempty(leadfield.leadfield{neighbors{j}(k)})
                chan_signal{i} = chan_signal{i} + leadfield.leadfield{neighbors{j}(k)}*current_mom{j}; %% just sum of signal. Should I use L2-norm?
            end
        end
    end
end
data_sim.trial = chan_signal; %% Simulated channel data.
if strcmp(data_sim.grad.type, 'bti248')
    datapow = [];
    tmpcfg = [];
    tmpcfg.output = 'pow';
    tmpcfg.method = 'mtmfft';
    tmpcfg.taper = 'dpss';
    tmpcfg.tapsmofrq =1;
    tmpcfg.keeptrials = 'no';
    datapow = ft_freqanalysis(tmpcfg, data_sim);
    
    subplot(2,3,3);
    tmpcfg = [];
    tmpcfg.layout = '4D248_helmet.mat';
    title('Simulated Channel Power');
    ft_topoplotER(tmpcfg, datapow);
end

subplot(2,3,4);
plot(data_sim.time{1}, data_sim.trial{1});
title('Simulated Channel Signals');
clear i j current_mom neighbors neighbors_1 source_mom source_ori source_sim
clear amplitude chan_signal k roi signal tmpcfg

%% Invert and Display
% Invert with LCMV
tmpcfg                   = [];
tmpcfg.method            = 'lcmv';
% mesh.vertices            = leadfield.pos;
% mesh.faces               = leadfield.tri;
% leadfield.mom            = spm_mesh_normals(mesh)';
tmpcfg.grid              = leadfield;
tmpcfg.headmodel         = headmodel;
tmpcfg.lcmv.lambda = '10%';
tmpcfg.lcmv.fixedori = 'no';
tmpcfg.lcmv.projectnoise = 'yes';
tmpcfg.lcmv.keepleadfield = 'yes';
source_lcmv = ft_sourceanalysis(tmpcfg, data_sim);
source_lcmv = ft_sourcedescriptives([], source_lcmv);

scale_factor = median(floor(log10(source_lcmv.avg.nai(find(~isnan(source_lcmv.avg.nai))))));
source_lcmv.avg.nai = 10.^(log10(source_lcmv.avg.nai)-scale_factor);

% plot the neural activity index (power/noise)
subplot(2,3,5);
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
tmpcfg.eloreta.fixedori = 'no';
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