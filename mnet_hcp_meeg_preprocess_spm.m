function [D,fieldtripData] = mnet_hcp_meeg_preprocess_spm(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEEG prprocessing for HCP dataset by using SPM.                         %
%                                                                         %
% High-pass filter -> Downsample -> Low-pass filter -> Artefect Detection %
% Remove all intermediate spm files                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.05.14 12:12 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end
%% Set Path
rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
%%
try
    D = spm_eeg_load(fullfile(config.savePath,['affdspm_', rMEGRawData]));
    fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
catch
    fieldtripData = load(fullfile(rMEGPath, rMEGRawData));
    try
        D = spm_eeg_load(config.savePath, ['spm_',rMEGRawData]);
    catch
        D = spm_eeg_ft2spm(fieldtripData.data, fullfile(config.savePath,['spm_' ,rMEGRawData]));
    end
    %% Bandpass filter & downsampling
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = config.downsample;
    D = spm_eeg_downsample(S); % resample of time has some edge effect problem
    
    cfg = [];
    cfg.resamplefs = config.downsample;
    cfg.resamplemethod = 'resample';
    fieldtripData.data = ft_resampledata(cfg, fieldtripData.data);
    
    S = [];
    S.D = D;
    S.band = 'high';
    S.freq = config.freqband(1);
    D = spm_eeg_filter(S);
    
    S = [];
    S.D = D;
    S.band = 'low';
    S.freq = config.freqband(2);
    D = spm_eeg_filter(S);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = config.freqband;
    fieldtripData.data = ft_preprocessing(cfg, fieldtripData.data);
    %% Artefact detaction
    S = [];
    S.D = D;
    S.mode = 'reject';
    S.methods.channels = {'MEG'};
    S.methods.fun = 'threshchan';
    S.methods.settings.threshold = 3000;%(3000 fT)
    S.methods.settings.excwin = 1000;
    D = spm_eeg_artefact(S);
    delete(fullfile(config.savePath,['ffdspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['ffdspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['fdspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['fdspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['dspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['dspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['spm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['spm_', rMEGRawData, '.dat']));
end