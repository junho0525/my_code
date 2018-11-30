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
%% invert using fieldtrip Beamformer LCMV
%% Simulation with wavelet
% simulation positions
use_my_library('spm',0);
use_my_library('ft',1);

cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel2d;
cfg.grad = data.grad;
cfg.channel = data.label;

leadfield = ft_prepare_leadfield(cfg);

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, data);

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
source_lcmv = ft_sourceanalysis(cfg, data);
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
source_eloreta = ft_sourceanalysis(cfg, data);
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

%% Parcellate
%  load MNET AAL MAP
use_my_library('ft',0);
use_my_library('spm',1);
cd('~/matlabwork/mnet_new/mnet0.92/labels/');
load('monet_aal128.mat');
aalmap = spm_vol('monet_aal128.nii');
dim = aalmap.dim;
toMNI = aalmap.mat;
V = spm_read_vols(aalmap);
V = sparse(V(:)); % Sparse vectorize for memory efficiency
voxels = find(V);
voxels_coord = zeros(length(voxels),4);
voxels_coord(:,4) = 1; 
for i = 1:length(voxels)
    voxels_coord(i,3) = ceil(voxels(i)/(dim(1)*dim(2)));% z
    if ~mod(voxels(i),dim(1)*dim(2))
        voxels_coord(i,2) = ceil(dim(1)*dim(2)/dim(1));
    else
        voxels_coord(i,2) = ceil(mod(voxels(i),dim(1)*dim(2))/dim(1));
    end
    if ~mod(voxels(i),dim(1))
        voxels_coord(i,1) = dim(1);
    else
        voxels_coord(i,1) = mod(voxels(i),dim(1));
    end
end
voxels_coord = toMNI*voxels_coord';
voxels_coord = voxels_coord';
voxels_coord(:,4) = [];
aalmap = [];
aalmap.label = AMAP.LABEL;
aalmap.pnt = voxels_coord;
aalmap.map = full(V(voxels));
aalmap.dim = dim;
cortex = find(aalmap.map<=90); % only use cortex
aalmap.map = aalmap.map(cortex);
aalmap.pnt = aalmap.pnt(cortex,:);
aalmap.label(91:end) = [];
clear AMAP voxels V voxels_coord dim i toMNI scale_factor cortex

cd(workingDir);