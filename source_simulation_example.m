%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG Source Simulation Example Script                                    %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.09.05 12:02 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load
addpath('/home/kuro/Codes/M_code/my_code');
%cd(workingDir);
use_my_library('ft',0);
use_my_library('spm',1);
load ('/projects1/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat');
headmodelfile = '/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_headmodel';
sourcemodelfile = '/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_sourcemodel_2d';
transformfile = '/projects1/HCPMEG/100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy/100307_MEG_anatomy_transform.txt';
%% Downsample data
if data.fsample>250
    cfg = [];
    cfg.resamplefs = 200;
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    data = ft_resampledata(cfg, data);
end
%% Simulation
cfg = [];
cfg.time = 2; % length of simulation
cfg.fsample = 200; % sampling frequency
cfg.wave_start = 0.5; % The time that sine wave starts.
cfg.wave_freq = 4; % Sine wave frequency
cfg.wave_amp  = 10; % Sime wave amplitude
cfg.wavelet_length = 0.5; % wavelet wavelength
cfg.wavelet_repetition = 1; % repetition of wavelet
cfg.noise_amp = 2; % noise amplitude

cfg.ROI = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIname = {'lLP', 'rLP', 'Prec', 'mPFC'};

cfg.ntrial = 1;
cfg.ROI_ori = {}; % Default : current randomly orients neighbor vertex.

mnet_source_simulation(cfg, data, sourcemodelfile, headmodelfile, transformfile);

%% Extract Source Timeseries
cfg = [];
cfg.method = 'lcmv';
cfg.resample_trials = 'rawtrial'; % 'no' - average all trials(Default), 'rawtrial' - do not do any resampling

source = mnet_source_reconstruction(cfg, data, sourcemodelfile, headmodelfile);
