addpath('/projects4/pi/jhson/my_code');
use_my_library('ft',0);
use_my_library('spm',1);

cfg.workingDir                          = '/projects4/pi/jhson/MEGResult';
cfg.savePath                            = '/projects4/pi/jhson/MEGResult/Restin';

% example files
headmodelfile = '/projects4/pi/jhson/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_headmodel.mat';
sourcemodelfile = '/projects4/pi/jhson/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_sourcemodel_2d.mat';
transformfile = '/projects4/pi/jhson/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_transform.txt';
fieldtripData = load('/projects4/pi/jhson/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc');
D =  spm_eeg_ft2spm(fieldtripData.data, fullfile(cfg.savePath,'spm_100307_MEG_3-Restin_rmegpreproc'));

cfg.preprocess.freqband                 = [0.1,40];
cfg.preprocess.downsample               = 200;
cfg.definetrial.bandfreq                = [8,12];

cfg.DCM.TimeWin                  = [0 2000];
cfg.DCM.FreqWin                  = [1,40];
cfg.DCM.analysis                 = 'CSD';
cfg.DCM.model                    = 'ERP';
cfg.DCM.spatial                  = 'ECD';
cfg.DCM.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.DCM.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.DCM.chanMode                 = 8;
cfg.DCM.trials                   = {[1], [2]}; % 1: Low alpha, 2: High Alpha
% Full Model
cfg.DCM.A{1}                     = ones(numel(cfg.DCM.ROIName),numel(cfg.DCM.ROIName)); % Forward
cfg.DCM.A{2}                     = ones(numel(cfg.DCM.ROIName),numel(cfg.DCM.ROIName)); % Backward
cfg.DCM.A{3}                     = ones(numel(cfg.DCM.ROIName),numel(cfg.DCM.ROIName)); % Modulatory

cfg.DCM.B{1}                     = zeros(numel(cfg.DCM.ROIName),numel(cfg.DCM.ROIName));

cfg.DCM.effect                   = zeros(1,0); % Between trial effect (Design matrix)
cfg.DCM.effectName               = {};

if ~isdir(fullfile(cfg.savePath,'DCMs'))
    mkdir(fullfile(cfg.savePath,'DCMs'));
end


[D, DCM_low, DCM_high] = mnet_HCPrMEG_AlphaContrastDCM(cfg, D, fieldtripData, sourcemodelfile, headmodelfile,transformfile);