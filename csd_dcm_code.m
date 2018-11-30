%% DCM %%
%% Specify parameters
%% Prepare data ( Set DCM.xY )
DCM = [];
DCM.xY.Dfile = D.fullfile;
DCM.xY.modality = 'MEG';
DCM.name = ['DCM_' D.fname];
DCM.val =D.val;
DCM.options.analysis = 'CSD';
DCM.options.model = 'ERP';
DCM.options.spatial = 'ECD';
DCM.options.trials = [1];
DCM.options.Tdcm = [1000 2000];
DCM.options.Fdcm = [1 100];
DCM.options.Nmodes = 8;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
DCM = spm_dcm_erp_data(DCM, 0);
spm_dcm_erp_results(DCM, 'Data');

%% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)
DCM.Lpos   = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
DCM.Sname  = {'lLP', 'rLP', 'Prec', 'mPFC'};
Nareas = size(DCM.Lpos,2);
% Spatial Model
DCM = spm_dcm_erp_dipfit_jhs(DCM);

%% Prepare Connectivity Matrices ( Set A, B, C - No need C(External input) for 'CSD' )
DCM.xU.name = {};
DCM.xU.X    = zeros(1,0);
% Specify connectivity model ( Set as full model )

DCM.A{1} = ones(Nareas, Nareas); % A{1} : Forward
DCM.A{2} = ones(Nareas,Nareas); % A{2} : Backward
DCM.A{3} = zeros(Nareas,Nareas); % A{3} : Modulatory
DCM.B{1} = zeros(Nareas,Nareas); % B{1} : Effect dependent modulatory
%% Invert
DCM = spm_dcm_csd_jhs(DCM);
