%% Set HCP MEG Resting Data Path
rawfile = '/media/kuro/DATA1/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat';
[rmegpath rawfile]    = fileparts(rawfile);
subjectid = strtok(rawfile, '_');
strind = strfind(rmegpath, subjectid);
megrootpath = rmegpath(1:(strind(1)-2)); % '/media/kuro/DATA1/HCPMEG'
rmegpath = rmegpath(strind(1):end); % '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc'
anapath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_anatomy'],subjectid, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
headmodelfile = [subjectid '_MEG_anatomy_headmodel.mat'];
sourcemodelfile = [subjectid '_MEG_anatomy_sourcemodel_2d.mat'];
clear strind;
%% Load Files
load(fullfile(megrootpath, rmegpath, rawfile)); % load HCP rMEG data
layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
headmodel = load(fullfile(megrootpath, anapath, headmodelfile));
headmodel = headmodel.headmodel;
sourcemodel = load(fullfile(megrootpath, anapath ,sourcemodelfile));
sourcemodel = sourcemodel.sourcemodel2d;
clear sourcemodelfile headmodelfile rawfile;

%% Calculate leadfield Without Data % channel number is not equal to data.label
cfg                        = [];
cfg.grad                   = data.grad; % sensor position
cfg.channel                = {'MEG'};
cfg.grid                   = sourcemodel; % sourcemodel
cfg.headmodel              = headmodel; % volume conduction model
lfwodata                   = ft_prepare_leadfield(cfg); % We can calculate leadfield without data itself.

%% Calculate Leadfield With Data
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = headmodel;
cfg.channel = {'MEG'};
lfwdata = ft_prepare_leadfield(cfg, data);
