%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HMM-MAR for HCP Resting State MEG                                       %
%                                                                         %
% This code is based on https://github.com/OHBA-analysis/HMM-MAR/blob/master/examples/NatComms2018_fullpipeline.m         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.06.26 16:17 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
addpath('/home/kuro/Codes/M_code/my_code');
addpath(genpath('/home/kuro/matlabwork/OHBA'));
mPath='/projects1/HCPMEG';
megDirInfo = dir(mPath);
subject = [];
for i=1:length(megDirInfo)
    if ~isempty(strfind(megDirInfo(i).name, 'MEG')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(megDirInfo(i).name,'_');
        subject(end).dir = megDirInfo(i).name;
    elseif ~isempty(str2num(megDirInfo(i).name)) % if subject folder is '[HCP data path]/[subject id]'
        subject(end+1).id = megDirInfo(i).name;
        subject(end).dir = megDirInfo(i).name;
    end
end
clear megDirInfo;
for i =1:length(subject)
    subjDirInfo = dir(fullfile(mPath,subject(i).dir));
    subject(i).anatomy = [];
    subject(i).Restin = [];
    subject(i).Wrkmem = [];
    subject(i).StoryM = [];
    subject(i).Motort = [];
    for j =1:length(subjDirInfo)
       if ~isempty(strfind(subjDirInfo(j).name, 'anatomy'))
           subject(i).anatomy = 1; 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Restin_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Restin/rmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Restin(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Restin =unique(subject(i).Restin); 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Wrkmem_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Wrkmem/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Wrkmem(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Wrkmem =unique(subject(i).Wrkmem);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'StoryM_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/StoryM/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).StoryM(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).StoryM =unique(subject(i).StoryM);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Motort_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Motort/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Motort(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Motort =unique(subject(i).Motort);
       end
    end
end
clear tempDirInfo subjDirInfo subDirInfo sessionNum i j k;

workingDir = '/projects1/MEGResult/HMM_MAR';

%% load
cd(workingDir);
use_my_library('ft',0);
use_my_library('spm',1);
load ('/projects1/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat');
load ('/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_headmodel');
load ('/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_sourcemodel_2d');

% data.trial(2:end)= [];
% data.time(2:end)= [];
% data.trialinfo(2:end,:)=[];
% D = spm_eeg_ft2spm(fieldtripData, fullfile(workingDir, 'spm_HMM_MAR_test'));
% D = spm_eeg_load('spm_100307_3_HMM_MAR');
%% Downsample data
if data.fsample>250
    cfg = [];
    cfg.resamplefs = 200;
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    data = ft_resampledata(cfg, data);
end
%% Source localize using OHBA Beamformer LCMV
cfg = [];
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodel2d;

[source ,weights, leadfield, cfg]= bf_inverse_lcmv_jhs(data,cfg);

% Plot Source on Mesh
for i = 1:size(source, 3)
    source_pow.pow(:,i) = sum(abs(fft(source(:,:,i)')),1)';
end
source_pow.pow = sum(source_pow.pow,2);
source_pow.pos = sourcemodel2d.pos;
source_pow.tri = sourcemodel2d.tri;
source_pow.inside = leadfield.inside;
source_pow.mask = (1+tanh(2.*(source_pow.pow./max(source_pow.pow(:))-0.5)))./2; 

% Source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_pow);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

%% Simulation
% Show selected source position
target = 3;
neighbors = find_neighbors_jhs(target, sourcemodel2d.tri,2);

bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;

source_sim = zeros(8004,1);
source_sim(neighbors) = 10;

ft_plot_mesh(bnd, 'vertexcolor', source_sim);
light;
% Simulation Data, Source Localize Using LCMV, and Source Plot
cfg = [];
cfg.vol = headmodel;
cfg.grad = data.grad;
cfg.fsample = data.fsample;
cfg.triallength = 2;
cfg.dip.pos = leadfield.pos(target,:);    % you can vary the location, here the dipole is along the z-axis
cfg.dip.mom = [1 0 0]';   % the dipole points along the x-axis
cfg.relnoise = 0.1;
cfg.ntrials = 147;
cfg.channel = data.label;
y0 = ft_dipolesimulation(cfg);

cfg = [];
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodel2d;

[source_sim ,weights_sim, leadfield_sim]= bf_inverse_lcmv_jhs(y0,cfg);

for i = 1:size(source_sim, 3)
    source_pow.pow(:,i) = sum(abs(fft(source_sim(:,:,i)')),1)';
end
source_pow.pow = sum(source_pow.pow,2);
source_pow.pos = sourcemodel2d.pos;
source_pow.tri = sourcemodel2d.tri;
source_pow.inside = leadfield_sim.inside;
sou rce_pow.mask = (1+tanh(2.*(source_pow.pow./max(source_pow.pow(:))-0.5)))./2; 

% Source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_pow);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

%% Source Localizing using fieldtrip beamformer lcmv
use_my_library('spm',0);
use_my_library('ft',1);
lcmv.lambda = '10%';
lcmv.fixedori = 'yes';
lcmv.porjectnoise = 'yes';
lcmv.keepleadfield = 'yes';
lcmvparam = ft_cfg2keyval(lcmv);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

data.leadfield = ft_prepare_leadfield(cfg);

source = beamformer_lcmv_jhs(sourcemodel2d, data.grad, headmodel, data, eye(length(data.label),length(data.label)), lcmvparam{:})
% Plot Source on Mesh
source_pow.pow = source{1}.pow;
source_pow.pos = sourcemodel2d.pos;
source_pow.tri = sourcemodel2d.tri;
source_pow.inside = source{1}.inside;
source_pow.mask = (1+tanh(2.*(source_pow.pow./max(source_pow.pow(:))-0.5)))./2; 

% Source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_pow);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

%% Simulation of fieldtrip beamformer lcmv
% Use same source position used above.
% Simulation Data, Source Localize Using LCMV, and Source Plot
use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = y0.grad;
cfg.channel = y0.label;

y0.leadfield = ft_prepare_leadfield(cfg);

source_sim = beamformer_lcmv_jhs(sourcemodel2d, y0.grad, headmodel, y0, eye(length(y0.label),length(y0.label)), lcmvparam{:})

% [source_sim ,weights_sim, leadfield_sim]= bf_inverse_lcmv_jhs(y0,cfg);
% Plot Source on Mesh
source_pow = [];
source_pow.pow = source_sim{1}.pow;
source_pow.pos = sourcemodel2d.pos;
source_pow.tri = sourcemodel2d.tri;
source_pow.inside = source_sim{1}.inside;
source_pow.mask = (1+tanh(2.*(source_pow.pow./max(source_pow.pow(:))-0.5)))./2; 

% Source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_pow);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

source_pow_des = ft_sourcedescriptives([], source_pow); % to get the neural-activity-index
% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0 8]; % [0.0 8]
cfg.opacitylim    = [1 8]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_pow_des);
view([-90 30]);
light;


%% Source Localizeing using frequency domain source localization methods(pcc)
%% Frequency Analysis
cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.fio = 10;
%cfg.foilim = [1 30]; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, data);

cfg                   = [];
cfg.frequency         = 10;
cfg.method            = 'dics';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.dics.lambda       = '10%';
cfg.dics.projectnoise = 'yes';
cfg.dics.keepfilter   = 'yes';
source = ft_sourceanalysis(cfg, spectdata);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [1 8]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);

%%
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepleadfield = 'yes';
cfg.singletrial = 'yes'
source_lcmv = ft_sourceanalysis(cfg, data);
source_lcmv = ft_sourcedescriptives([], source_lcmv);

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';

source_lcmv_rescale = source_lcmv;
source_lcmv_rescale.avg.nai = 10.^(log10(source_lcmv_rescale.avg.nai)-floor(log10(source_lcmv_rescale.avg.nai)));

ft_sourceplot(cfg, source_lcmv_rescale);

% Show selected source position
target = 3;
neighbors = find_neighbors_jhs(target, sourcemodel2d.tri,2);

bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;

source_sim = zeros(8004,1);
source_sim(neighbors) = 10;

ft_plot_mesh(bnd, 'vertexcolor', source_sim);
light;
% Simulation Data, Source Localize Using LCMV, and Source Plot
cfg = [];
cfg.vol = headmodel;
cfg.grad = data.grad;
cfg.fsample = data.fsample;
cfg.triallength = 2;
cfg.dip.pos = leadfield.pos(target,:);    % you can vary the location, here the dipole is along the z-axis
cfg.dip.mom = [1 0 0]';   % the dipole points along the x-axis
cfg.dip.amplitude = [10]
cfg.relnoise = 0.1;
cfg.ntrials = 147;
cfg.channel = data.label;
y0 = ft_dipolesimulation(cfg);

cfg                   = [];
cfg.method            = 'lcmv';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepleadfield = 'yes';
cfg.singletrial = 'yes'
source_lcmv_sim = ft_sourceanalysis(cfg, y0);
source_lcmv_sim = ft_sourcedescriptives([], source_lcmv_sim);

source_lcmv_sim = source_lcmv_sim;
source_lcmv_sim.avg.nai = 10.^(log10(source_lcmv_sim.avg.nai)-floor(log10(source_lcmv_sim.avg.nai)));

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_lcmv_sim);
%% Simulation
t = data.time{1};
signal = 10*sin(2*pi*10*t)+(rand(size(t))-0.5)*4;
target = 300;
neighbors = find_neighbors_jhs(target, sourcemodel2d.tri,2);

bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;

source_sim = zeros(8004,1);
source_sim(neighbors) = sum(abs(fft(signal)));

ft_plot_mesh(bnd, 'vertexcolor', source_sim);
light;

neighbors_1 = setdiff(find_neighbors_jhs(target, sourcemodel2d.tri,1),[3]);
use_my_library('spm',0);
use_my_library('ft',1);
source_ori =(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(3,:))/norm(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(3,:));
source_mom = [];
sensor_signal = [];
for i=1:147
    signal = 10*sin(2*pi*10*t)+(rand(size(t))-0.5)*4;
    current_mom(1,:) = source_ori(1).*signal;
    current_mom(2,:) = source_ori(2).*signal;
    current_mom(3,:) = source_ori(3).*signal;
    source_mom{i} = current_mom;
    chan_signal{i} = zeros(241,401);
    for j = 1:length(neighbors)
        chan_signal{i} = chan_signal{i} + leadfield.leadfield{neighbors(j)}*source_mom{i}; %% just sum of signal. Should I use L2-norm?
    end
end

data_sim = data;
data_sim.trial = chan_signal; %% Simulated channel data.

cfg                   = [];
cfg.method            = 'lcmv';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepleadfield = 'yes';
cfg.singletrial = 'yes'
source_lcmv_sim = ft_sourceanalysis(cfg, data_sim);
source_lcmv_sim = ft_sourcedescriptives([], source_lcmv_sim);

source_lcmv_sim = source_lcmv_sim;
source_lcmv_sim.avg.nai = 10.^(log10(source_lcmv_sim.avg.nai)-floor(log10(source_lcmv_sim.avg.nai)));

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_lcmv_sim);

%% Simulation with SPM
use_my_library('ft',0);
use_my_library('spm',1);
D = spm_eeg_ft2spm(data, 'spm_HMM_MAR_test'); 

%% Specify Forward Model
% read the source and volume conduction model from current dir with
% outputs of previous pipelines
fid = fopen(fullfile('/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_transform.txt'), 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
fid = [];
clear fid;
strFid = [];
clear strFid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify sourcemodel, headmodel, sensor and transform if it necessary
sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
sourcemodelsubj = sourcemodel2d;
headmodel = ft_convert_units(headmodel, 'cm');
grad=data.grad;
grad = ft_convert_units(grad, 'cm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the component data in order for ft_sourceanalysis to be able to
% swallow it
channels = data.label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the forward solution

job = [];
job.vol = headmodel;
job.grid = sourcemodelsubj;
job.grad = grad;
job.channel = channels;
job.normalize = 'yes';
job.reducerank = 2;

val = []; tlck=[]; sourcemodelsubj = []; resultprefix = []; options = []; mixing = []; megpath = []; rmeg = []; mainpath = [];
i = []; grad = []; ft_default = []; experimentid = []; comp_class = []; channels = []; anaPath = [];
clear val tlck sourcemodelsubj resultprefix options mixing megpath;
clear rmeg mainpath i grad ft_default experimentid comp_class channels anaPath;

% specify forward model (Bring forward model information 
D.inv = [];
D.inv{1}.forward = [];
D.inv{1}.method='Imaging';
D.inv{1}.forward.toMNI = transform.bti2spm;
D.inv{1}.forward.fromMNI = transform.spm2bti;
D.inv{1}.forward.voltype='Single Shell';
D.inv{1}.forward.vol = ft_convert_units(job.vol,'m');
D.inv{1}.forward.modality = 'MEG';
D.inv{1}.forward.siunits=1;
D.inv{1}.forward.mesh=[];
D.inv{1}.forward.mesh_correction=[];
D.inv{1}.forward.mesh.face = job.grid.tri;
grid = ft_convert_units(job.grid,'m');
D.inv{1}.forward.mesh.vert = grid.pos;
D.inv{1}.forward.sensors = [];
D.inv{1}.forward.sensors = ft_convert_units(job.grad,'m');
D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face;
grid = ft_convert_units(job.grid,'mm');
D.inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(D.inv{1}.forward.toMNI,grid.pos); % tess_mni should have mni coordinate
D.inv{1}.forward.toMNI = [D.inv{1}.forward.toMNI(:, 1:3)*1000 D.inv{1}.forward.toMNI(:,4)];
D.inv{1}.forward.fromMNI = eye(4)/D.inv{1}.forward.toMNI
D.inv{1}.datareg.modality = 'MEG';
D.inv{1}.datareg.sensors = ft_convert_units(job.grad,'m');
D.inv{1}.datareg.toMNI = D.inv{1}.forward.toMNI;
D.inv{1}.datareg.fromMNI = D.inv{1}.forward.fromMNI;
Datareg = []; grid = [];
save(fullfile(D.path,D.fname),'D');
clear Datareg grid job;
% check meshes and co-registration
% if config.showForwardResult; spm_eeg_inv_checkforward(D); end

%% Simulation Code
% specify job for spm_cfg_eeg_inv_simulate.m, function run_simulation(job)
use_my_library('spm',0);
use_my_library('ft',1);

t = data.time{1};
signal = 10*sin(2*pi*10*t)+(rand(size(t))-0.5)*4;

target = [1608 862 2839 2078 3587 7181 7658 6110];
figure;
suptitle('suptitle')
for i = 1:length(target)
    subplot(2,4,i);
    neighbors = find_neighbors_jhs(target(i), sourcemodel2d.tri,2);

    bnd.pnt = sourcemodel2d.pos;
    bnd.tri = sourcemodel2d.tri;

    source_sim = zeros(8004,1);
    source_sim(neighbors) = sum(abs(fft(signal)));

    ft_plot_mesh(bnd, 'vertexcolor', source_sim);
    light;
    title(num2str(target(i)));
end

%%
use_my_library('ft',0);
use_my_library('spm',1);

job = [];
job.D = [];
fname = strsplit(D.fname, '.');
job.D{1} = fname{1};
job.val = 1;
job.prefix = 'sim_';
job.whatconditions.all =1;
job.isinversion.setsources.woi = [0 2000]; % Specify Time Window that source simulation starts and ends
job.isinversion.setsources.dipmom=[10 5]; % Specify Dipole Moment
% Three options
% dipmom = [intensity] or [intensity, smoothness_radius_in_mm], [x, y, z]
job.isinversion.setsources.isSin.foi = 10;
job.isinversion.setsources.locs = D.inv{1}.mesh.tess_mni.vert(target(1),:);
% [-6 28 34; -6 -52 40; -46 -70 36 ;-48 -69 7; -39 -7 56; -28 52 19]; %ACC, PCC, LLP,
% LV1, LMT, LFEF, L aPFC
job.isSNR.setSNR = 1;

%
% Run Simulation
% 
% modified code of spm_cfg_eeg_inv_simulate.m
%
trialind=[];
D.con = 1;
if isfield(job.isSNR,'whitenoise'),
    whitenoisefT=job.isSNR.whitenoise; %% internal units as femto Tesla
    SNRdB=[];
else
    SNRdB=job.isSNR.setSNR;
    whitenoisefT=[];
end;

Nsources=size(job.isinversion.setsources.locs,1)
mnimesh=[]; %% using mesh defined in forward model at the moment
ormni=[]; %% dipoles will get orientations from cortical mesh
dipmom=job.isinversion.setsources.dipmom;
if size(job.isinversion.setsources.dipmom,2)==3, %% dipole moment has 3 dimension
    disp('Simulating dipoles without reference to cortical surface');
    for j=1:size(dipmom,1),
        ormni(j,:)=dipmom(j,:)./sqrt(dot(dipmom(j,:),dipmom(j,:))); %% unit orientation in MNI space
        nAmdipmom(j)=sqrt(dot(dipmom(j,:),dipmom(j,:))); % magnitude of dipole
    end;
    dipfwhm = [];
else %% only one moment parameter given
    nAmdipmom=dipmom(:,1); %% total momnent in nAm
    dipfwhm=[];
    if size(dipmom,1)==2,
        dipfwhm=dipmom(:,2); %% fhwm in mm
    end;
end;

woi=job.isinversion.setsources.woi./1000;
timeind = intersect(find(D.time>woi(1)),find(D.time<=woi(2)));
simsignal = zeros(Nsources,length(timeind));

if isfield(job.isinversion.setsources.isSin,'fband'),
    % Simulate orthogonal Gaussian signals

    simsignal=randn(Nsources,length(timeind));
    % filter to bandwidth
    simsignal=ft_preproc_lowpassfilter(simsignal,D.fsample,job.isinversion.setsources.isSin.fband(2),2);
    simsignal=ft_preproc_highpassfilter(simsignal,D.fsample,job.isinversion.setsources.isSin.fband(1),2);
    [u,s,v]=svd(simsignal*simsignal');
    simsignal=u*simsignal; %% orthogonalise all signals
end; % if isfield fband

%  simsignal=simsignal.*repmat(nAmdipmom,1,size(simsignal,2)); %% now scale by moment

if isfield(job.isinversion.setsources.isSin,'foi'),
    % simulate sinusoids
    sin_freq=job.isinversion.setsources.isSin.foi;
    % Create the waveform for each source

    for j=1:Nsources                % For each source
        simsignal(j,:)=sin((D.time(timeind)- D.time(min(timeind)))*sin_freq(j)*2*pi);
    end; % for j


end; %% if isfield foi

simsignal=simsignal./repmat(std(simsignal'),size(simsignal,2),1)'; %% Set sim signal to have unit variance
%[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.isinversion.setsources.locs,simsignal,woi,whitenoisefT,SNRdB,trialind,mnimesh,SmthInit);
figure();
plot(D.time(timeind),simsignal);
xlabel('time');
ylabel('Normalized amplitude (to be scaled by dip moment later)');
legend(num2str([1:Nsources]'));
[simD,meshsourceind]=spm_eeg_simulate({D},job.prefix,job.isinversion.setsources.locs,simsignal,ormni,woi,whitenoisefT,SNRdB,trialind,mnimesh,dipfwhm,nAmdipmom);
save(D);
clear job trialind fname whitenoisefT SNRdB Nsources mnimesh ormni dipmom nAmdipmom dipfwhm woi timeind sinfreq
%% Invert using fieldtrip mne and plot using spm_mip
% invert using fieldtrip mne
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
cfg.demean = 'yes';
sim_data = ft_resampledata(cfg, data);

trial = {};
for i = 1:size(simD(:,:,:),3)
    trial{i} = simD(:,:,i);
end
sim_data.trial = trial;
clear trial;

use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg                   = [];
cfg.method            = 'mne';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.mne.lambda = '10%';
cfg.mne.fixedori = 'yes';
cfg.mne.projectnoise = 'yes';
cfg.mne.keepleadfield = 'yes';
%cfg.singletrial = 'yes'
source = ft_sourceanalysis(cfg, sim_data);
source = ft_sourcedescriptives([], source);

use_my_library('ft',0);
use_my_library('spm',1);
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
[T,i] = sort(-abs(source.avg.pow));
i = i(1:512);
spm_mip(source.avg.pow(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image

source.avg.pow = sum(moment,2);

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylom    = [0.2 0.5];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);


%% Invert using fieldtrip eloreta and plot using spm_mip
% invert using fieldtrip eloreta
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
cfg.demean = 'yes';
sim_data = ft_resampledata(cfg, data);

trial = {};
for i = 1:size(simD(:,:,:),3)
    trial{i} = simD(:,:,i);
end
sim_data.trial = trial;
clear trial;

use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg                   = [];
cfg.method            = 'eloreta';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.eloreta.lambda = 10;
cfg.eloreta.fixedori = 'yes';
cfg.eloreta.projectnoise = 'yes';
cfg.eloreta.keepleadfield = 'yes';
%cfg.singletrial = 'yes'
source = ft_sourceanalysis(cfg, sim_data);
source = ft_sourcedescriptives([], source);


% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);

use_my_library('ft',0);
use_my_library('spm',1);
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
moment = zeros(length(source.avg.mom),length(source.avg.mom{1}));
for i =  1:size(moment,1)
    if ~isempty(source.avg.mom{i})
        moment(i,:) = sum(abs(source.avg.mom{i}),1);
    end
end
[i,js] = max(max(abs(moment),[],2));
[i,jt]=max(abs(moment(js,:)));
Js = moment(:,jt);
[T,i] = sort(-abs(Js));
i = i(1:64);
spm_mip(source.avg.pow(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image

%% Invert using fieldtrip pcc and plot using spm_mip
% Frequency Analysis
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
cfg.demean = 'yes';
sim_data = ft_resampledata(cfg, data);

trial = {};
for i = 1:size(simD(:,:,:),3)
    trial{i} = simD(:,:,i);
end
sim_data.trial = trial;
clear trial;

use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foilim = [1,50];
%cfg.foilim = [1 30]; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, sim_data);

cfg                   = [];
cfg.frequency         = [1,50];
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.pcc.lambda       = '10%';
cfg.pcc.projectnoise = 'yes';
cfg.pcc.keepfilter   = 'yes';
source = ft_sourceanalysis(cfg, spectdata);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [1.5 1.7]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);

% plot the power
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = 'auto';
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

% plot nai using spm_mip
use_my_library('ft',0);
use_my_library('spm',1);
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
subplot(2,1,1);
[T,i] = sort(-abs(source.avg.nai));
i = i(1:512);
spm_mip(source.avg.nai(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image
title('nai')
subplot(2,1,2);
moment = zeros(length(source.avg.mom),length(source.avg.mom{1}));
for i =  1:size(moment,1)
    if ~isempty(source.avg.mom{i})
        moment(i,:) = sum(abs(source.avg.mom{i}),1);
    end
end
[i,js] = max(max(abs(moment),[],2));
[i,jt]=max(abs(moment(js,:)));
Js = moment(:,jt);
[T,i] = sort(-abs(Js));
i = i(1:512);
spm_mip(Js(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image
title('moment');


%% Invert using fieldtrip beamformer lcmv and plot using spm_mip
% invert using fieldtrip Beamformer LCMV
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
cfg.demean = 'yes';
sim_data = ft_resampledata(cfg, data);

trial = {};
for i = 1:size(simD(:,:,:),3)
    trial{i} = simD(:,:,i);
end
sim_data.trial = trial;
clear trial;

use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg                   = [];
cfg.method            = 'lcmv';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepleadfield = 'yes';
%cfg.singletrial = 'yes'
source_lcmv_sim = ft_sourceanalysis(cfg, sim_data);
source_lcmv_sim = ft_sourcedescriptives([], source_lcmv_sim);

source_lcmv_sim.avg.nai = 10.^(log10(source_lcmv_sim.avg.nai)-floor(log10(source_lcmv_sim.avg.nai)));

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_lcmv_sim);

use_my_library('ft',0);
use_my_library('spm',1);
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
subplot(2,1,1);
[T,i] = sort(-abs(source_lcmv_sim.avg.nai));
i = i(1:64);
spm_mip(source_lcmv_sim.avg.nai(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image
subplot(2,1,2);
moment = zeros(length(source_lcmv_sim.avg.mom),length(source_lcmv_sim.avg.mom{1}));
for i =  1:size(moment,1)
    if ~isempty(source_lcmv_sim.avg.mom{i})
        moment(i,:) = source_lcmv_sim.avg.mom{i};
    end
end
[i,js] = max(max(abs(moment),[],2));
[i,jt]=max(abs(moment(js,:)));
Js = moment(:,jt);
[T,i] = sort(-abs(Js));
i = i(1:512);
spm_mip(Js(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image

%% Invert using fieldtrip dics and plot using spm_mip
% Frequency Analysis
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
cfg.demean = 'yes';
sim_data = ft_resampledata(cfg, data);

trial = {};
for i = 1:size(simD(:,:,:),3)
    trial{i} = simD(:,:,i);
end
sim_data.trial = trial;
clear trial;

use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foilim = [1,50];
%cfg.foilim = [1 30]; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, sim_data);

cfg                   = [];
cfg.frequency         = [1,50];
cfg.method            = 'dics';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.dics.lambda       = '10%';
cfg.dics.projectnoise = 'yes';
cfg.dics.keepfilter   = 'yes';
source = ft_sourceanalysis(cfg, spectdata);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [.8 1]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);

% plot the power
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = 'auto';
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

% plot nai using spm_mip
use_my_library('ft',0);
use_my_library('spm',1);
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
% subplot(2,1,1);
[T,i] = sort(-abs(source.avg.nai));
i = i(1:64);
spm_mip(source.avg.nai(i),D.inv{1}.mesh.tess_mni.vert(i,:)',6);
axis image
title('nai');
clear i
%% Invert SPM Example
% Specify Inverting Methods.
% Use 'LOR', 'GS', 'IID', 'EEB' for Methods 'COH', 'MSP(GS)', 'IID', 'EBB'
use_my_library('ft',0);
use_my_library('spm',1);
simD.inv{1}.inverse.type = 'GS';
simD.inv{1}.inverse.woi = [0,2000];
simD.inv{1}.inverse.Han = 1;
simD.inv{1}.inverse.lpf = 1;
simD.inv{1}.inverse.hpf = 40;
simD.inv{1}.inverse.modality = {'MEG'};

%    inverse.woi  = fix([max(min(job.isstandard.custom.woi), 1000*D.time(1)) min(max(job.isstandard.custom.woi), 1000*D.time(end))]);

simD = spm_eeg_invert(simD);
save(simD);

% use source of eloreta
source.avg.pow =  full(sum(abs(sparse(simD.inv{1}.inverse.J{1}*simD.inv{1}.inverse.T')),2));

use_my_library('spm',0);
use_my_library('ft',1);

% source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolorlim   = 'auto';
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

%% Simulation with wavelet
% simulation positions
use_my_library('spm',0);
use_my_library('ft',1);

t = data.time{1}; % time for sin
a = 0.5; % The time that wavelet starts
sin_freq = 10; % sinwave frequency
sig_amp = 10; % peak signal amplitude
wave_length = 0.5; % wavelet wavelength
wavelet_repetition = 1; % repetition of wavelet
noise_amp = 2; % noise amplitude


signal = zeros(size(t));
signal(find((t>a).* (t<a+wave_length*wavelet_repetition))) = ... 
    (0.5 - 0.5*cos(2*pi*(t(find((t>a).* (t<a+wave_length*wavelet_repetition)))-a)/wave_length)) ... 
    .*(sig_amp*sin(2*pi*sin_freq*(t(find((t>a).* (t<a+wave_length*wavelet_repetition)))-a)));
signal = signal +2*(rand(size(t))-0.5)*noise_amp;
figure;
plot(t,signal);

target = [1608 862 2839 2078 3587 7181 7658 6110];
target = target(4);
figure;
neighbors = find_neighbors_jhs(target, sourcemodel2d.tri,2);

bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;
source_sim = zeros(8004,1);
source_sim(neighbors) = sum(abs(fft(signal)));

ft_plot_mesh(bnd, 'vertexcolor', source_sim);
light;

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

neighbors_1 = setdiff(find_neighbors_jhs(target, sourcemodel2d.tri,1),[3]);
use_my_library('spm',0);
use_my_library('ft',1);
source_ori =(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(3,:))/norm(sourcemodel2d.pos(neighbors_1(1),:) - sourcemodel2d.pos(3,:));
source_mom = [];
sensor_signal = [];
for i=1:147
    current_mom(1,:) = source_ori(1).*signal;
    current_mom(2,:) = source_ori(2).*signal;
    current_mom(3,:) = source_ori(3).*signal;
    source_mom{i} = current_mom;
    chan_signal{i} = zeros(241,401);
    for j = 1:length(neighbors)
        chan_signal{i} = chan_signal{i} + leadfield.leadfield{neighbors(j)}*source_mom{i}; %% just sum of signal. Should I use L2-norm?
    end
end

data_sim = data;
data_sim.trial = chan_signal; %% Simulated channel data.

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, data_sim);

clear i j bnd cfg current_mom neighbors neighbors_1 noise_amp sig_amp sin_freq t source_mom source_ori source_sim sinfreq wave_length
clear wavelet_repetition a amplitude chan_signal

% Invert
% Invert with LCMV
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepleadfield = 'yes';
%cfg.singletrial = 'yes'
source_lcmv = ft_sourceanalysis(cfg, data_sim);
source_lcmv = ft_sourcedescriptives([], source_lcmv);

scale_factor = median(floor(log10(source_lcmv.avg.nai(find(~isnan(source_lcmv.avg.nai))))));
source_lcmv.avg.nai = 10.^(log10(source_lcmv.avg.nai)-scale_factor);

% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = 'nai';
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_lcmv);
title('LCMV nai');

cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_lcmv);
title('LCMV pow');
% Invert with loreta
cfg                   = [];
cfg.method            = 'eloreta';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.eloreta.lambda = 10;
cfg.eloreta.fixedori = 'yes';
cfg.eloreta.projectnoise = 'yes';
cfg.eloreta.keepleadfield = 'yes';
source_eloreta = ft_sourceanalysis(cfg, data_sim);
source_eloreta = ft_sourcedescriptives([], source_eloreta);


% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source_eloreta);
title('eLORETA pow');


