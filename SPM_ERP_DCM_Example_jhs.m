%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This pipeline code is for spm MEG task DCM example                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.04.02 20:22 - By Junho Son                                     %
%         - Make this code based on the code rMEG_pipeline_HCP_jhs_v2.0.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DCM = [];
DCM.xY.Dfile = D.fullfile;
DCM.options.analysis = 'ERP'; % 'ERP', 'CSD', 'TFM', 'IND', 'PHA', 'NFM'
DCM.options.model    = 'ERP'; % 'ERP', 'SEP', 'CMC', 'LFP', 'NMM', 'MFM', 'CMM', 'NMDA', 'CMM_NMDA', 'NFM'
DCM.options.spatial  = 'ECD'; % 'ECD', 'IMG', 'LFP', 'ITR'
DCM.options.trials   = [1 2 3]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 500;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 8;     % nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
DCM.options.onset    = 60;    % selection of onset (prior mean)
DCM.options.D        = 1;     % downsampling

% Data and spatial model
DCM  = spm_dcm_erp_data(DCM);

% Location priors for dipoles
DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};
Nareas    = size(DCM.Lpos,2);

% Spatial model
DCM = spm_dcm_erp_dipfit_jhs(DCM);

% Specify connectivity model
DCM.A{1} = zeros(Nareas,Nareas); % A{1} : forward
DCM.A{1} = zeros(Nareas, Nareas);
DCM.A{1}(3,1) = 1;
DCM.A{1}(4,2) = 1;
DCM.A{1}(5,4) = 1;

DCM.A{2} = zeros(Nareas,Nareas); % A{2} : backward
DCM.A{2}(1,3) = 1;
DCM.A{2}(2,4) = 1;
DCM.A{2}(4,5) = 1;

DCM.A{3} = zeros(Nareas,Nareas); % A{3} : lateral
DCM.A{3}(4,3) = 1;
DCM.A{3}(3,4) = 1;

DCM.B{1} = DCM.A{1} + DCM.A{2}; % B{1} : between trial effect(trial dependent)
DCM.B{1}(1,1) = 1;
DCM.B{1}(2,2) = 1;

DCM.B{2} = DCM.A{1} + DCM.A{2}; % B{1} : between trial effect(trial dependent)
DCM.B{2}(1,1) = 1;
DCM.B{2}(2,2) = 1;

DCM.C = [1; 1; 0; 0; 0]; % input

% Between trial effects
DCM.xU.X = [1 0; 0 1; 0 0]
DCM.xU.name = {'Zero-back', 'Two-back'};

% Invert
DCM.name = 'DCMexample';
DCM      = spm_dcm_erp_jhs(DCM);

