%% 
load 100307_MEG_3-Restin_rmegpreproc.mat; % load fieldtrip data struct
D = spm_eeg_ft2spm(data, 'test_data');
fieldtripData = D.ftraw; % make converted fieldtrip data from spm data structure
fieldtripData.fsample = D.fsample;
%% Resampling for 200Hz
if fieldtripData.fsample >250
    cfg = [];
    cfg.resamplefs = 250;
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    fieldtripData = ft_resampledata(cfg, fieldtripData);
end
%% Band-pass filter
cfg = [];
cfg.bpfilter ='yes';
cfg.bpfreq = [0.1, 100]; %% must satisfy this condition ( freqmax*2 < fsample )
fieldtripData = ft_preprocessing(cfg, fieldtripData);
%% Data concatenation
tmpdata = [];
tmptime = [];
data_concat = fieldtripData;
for i = 1:length(fieldtripData.time)
    starttime = (i-1)*(length(fieldtripData.time{i})-1)+1;
    endtime = i*(length(fieldtripData.time{i})-1);
    tmpdata(:,starttime:endtime) = fieldtripData.trial{i}(:,1:length(fieldtripData.time{i})-1);
    tmptime(:,starttime:endtime) = fieldtripData.time{i}(:,1:length(fieldtripData.time{i})-1)+2*(i-1);
end
data_concat.time = {tmptime};
data_concat.trial = {tmpdata};
data_concat = rmfield(data_concat,{'trialinfo'});
clear tmpdata tmptime starttime endtime
%% ICA
cfg = [];
cfg.method = 'runica';
cfg.runica.lrate = 1e-3; % Learning rate
cfg.runica.pca = 20; % Number of components
data_comp = ft_componentanalysis(cfg,data_concat);
cfg = [];
cfg.channel =1:10;
cfg.viewmode = 'component';
cfg.layout = '4D248_helmet.mat';
ft_databrowser(cfg, data_comp);
%% frequency analysis of components
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'pow';
cfg.tapsmofrq = 1;
%cfg.keeptrials = 'yes';
cfg.foi = 2:30;
datapow = ft_freqanalysis(cfg, data_comp); % Error in this code
% ft_freqanalysis reports error in 510th line.
% Stop at 463th line, excute code below
% >> cfg.pad = ceil(cfg.pad);
% This seem to floating point problem. 
figure;
for i = 1:20
subplot(10,2,i);
cfg.channel = datapow.label{i};
ft_singleplotER(cfg, datapow);
end
clear datapow
%% Load headmodel & sourcemodel
load('E:\matlabworks\fieldtrip\fieldtrip-20180918\template\headmodel\standard_singleshell.mat'); % headmodel
sourcemodel = ft_read_headshape('E:\matlabworks\fieldtrip\fieldtrip-20180918\template\sourcemodel\cortex_8196.surf.gii');% load 2d sourcemodel
headmodel = vol;
clear vol
%% Construct Leadfield
cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.grad = fieldtripData.grad;
cfg.channel = fieldtripData.label;
leadfield = ft_prepare_leadfield(cfg);
%% Generate Data for All Components.
tmp                 = data_concat.trial{1};
tmp(~isfinite(tmp)) = 0; % ensure that the scaling is a finite value
scale = norm((tmp*tmp')./size(tmp,2)); clear tmp;
scale = sqrt(scale);
if scale ~=0
    tmp = data_comp.trial{1} ./scale;
    for i = 1:size(data_comp.trial{1},1)
        data{i} = data_concat;
        data{i}.trial{1} = zeros(1,size(data_comp.topo,2));
        data{i}.trial{1}(i) = 1;
        data{i}.trial{1} = scale * data_comp.topo* diag(data{i}.trial{1})*tmp;
    end
    clear tmp
end
%% Source Reconstruction (using Beamformer LCMV)
cfg = [];
cfg.grid = leadfield;
cfg.headmodel = headmodel;
cfg.method = 'eloreta';
cfg.eloreta.lambda = 10;
cfg.eloreta.fixedori = 'yes';
cfg.eloreta.projectnoise = 'yes';
cfg.eloreta.keepleadfield = 'yes';
%% Source localize for each ICA components
for i =6:20
    source = ft_sourceanalysis(cfg, data{i});
    source = ft_sourcedescriptives([], source);

    %% Plot Sources
    figure;
    bnd.pnt = sourcemodel.pos;
    bnd.tri = sourcemodel.tri;
    ft_plot_mesh(bnd);
    ft_plot_mesh(bnd, 'vertexcolor', source.avg.pow, 'facealpha',source.avg.pow,'alphalim',[0 max(source.avg.pow)],'colormap','jet');
    lighting gouraud
    camlight
    title(sprintf('Source Reconstruction of Component %d\n(eLORETA)',i));
    drawnow
    clear source
end
%% Source simulation for check channel pattern of occipital source

cfg = [];
cfg.time = 2; % length of simulation
cfg.fsample = 200; % sampling frequency
cfg.wave_start = 0.5; % The time that sine wave starts.
cfg.wave_freq = 4; % Sine wave frequency
cfg.wave_amp  = 10; % Sime wave amplitude
cfg.wavelet_length = 0.5; % wavelet wavelength
cfg.wavelet_repetition = 1; % repetition of wavelet
cfg.noise_amp = 2; % noise amplitude

% cfg.ROI = [[-28;-96;-6] [-36; -86 ; -8] [28;-96;-6] [44;-84;-10] [-2;-80;34]];
% cfg.ROIname = {'lV1','lV5','rV1','rV5','V1'};
% cfg.ROI = [[10;64;18],[-10;64;18] [10;40;50],[-10;40;50]];
% cfg.ROIname = {'rFrontal Pole','lFrontal Pole' 'rSuperior Frontal Gyrus','lSuperior Frontal Gyrus'};
% cfg.ROI = [[-12;-16;74],[-32;-28;70] [-52 ;-4;54] [-63;-14;30]];
% cfg.ROIname = {'rPMC','rPMC','rPMC','rPMC'};

cfg.ROI = [[-38;16;50]];
cfg.ROIname = {'lBroca'};

dip.n_seeds = 1;
dip.n_dip = 1;
dip.Mtp = 1;
dip.j{1} = zeros(3,1);
dip.loc{1} = [-38;16;50];
spm_eeg_inv_ecd_DrawDip('init',dip);

cfg.ntrial = 1;
cfg.ROI_ori = {[0;1;0]}; % Default : current randomly orients neighbor vertex.
mesh.vertices = sourcemodel.pos;
mesh.faces = sourcemodel.tri;
norm = spm_mesh_normals(mesh, 1);
cfg.ROI_ori = mat2cell(norm', 3,ones(1,length(norm)));

mnet_source_simulation_v2(cfg, data, sourcemodel, headmodel, eye(4));

%% 
M.vertices = sourcemodel.pos;
M.faces = sourcemodel.tri
sourceOri = spm_mesh_normals(M)
bnd.pnt = M.vertices
bnd.tri = M.faces
bnd2 = bnd.pnt+sourceOri;
ft_plot_mesh(bnd);
camlight
lighting gouraud
alpha(0.3)
alpha(0.2)
alpha(0.05)
for i=1:8196
hold on
plot3([bnd.pnt(i,1);bnd2(i,1)],[bnd.pnt(i,2) ;bnd2(i,2)],[bnd.pnt(i,3);bnd2(i,3)])
end
