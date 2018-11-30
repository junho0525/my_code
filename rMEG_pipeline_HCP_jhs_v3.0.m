%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              This Is M/EEG Resting State DCM Pipeline Code              %
%                              By Using SPM                               %
%                                                                         %
% We need five inputs                                                     %
% You need to specify these parameters.( All of them are string )         %
%     fileType        - File type               ('EEG', 'MEG' or 'HCP')   %
%     rawFile         - M/EEG raw data filename (eg. Subject2Rest.ds)     %
%     subjectID       - Subject ID for HCP data (eg. '100307')            %
%     sourceModelFile - Source model filename   (eg. sourcemodel4k.mat)   %
%     headModelFile   - Head model filename     (eg. headmodel.mat)       %
%     layoutFile      - Channel layout file     (eg. 4D248_helmet.mat)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.02.17 18:23 - By Junho Son                                     %
%     2018.05.01       - By Junho Son                                     %
%               Add DCM                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SET DATA PATH    %%
%% Check Input Parameters
%  If you want to use example datapath, please excute this section twice.
if ~exist('fileType')
   warning('You should specify file type.');
   fileType = 'HCP';
end
if exist('rawFile')
    isRawFile=1;
else
    isRawFile =0;
end
if exist('subjectID')
    isSubjectID = 1;
else
    isSubjectID = 0;
end

% distract path from rawFile
if ~isRawFile&&~isSubjectID
    warning('There is no rawFile. This code use default dataset (HCP-subject 100307 resting data)');
    % rawFile = 'D:\HCPDATA\100307_MEG\100307_MEG_Restin_preproc\100307\MEG\Restin\rmegpreproc\100307_MEG_3-Restin_rmegpreproc';
    rawFile = '/media/kuro/DATA1/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat';
end
switch fileType
    case 'HCP'
        if isRawFile
            [rMEGPath rawFile]    = fileparts(rawFile);
            subjectID = strtok(rawFile, '_');
            strind = strfind(rMEGPath, subjectID);
            
            megRootPath = rMEGPath(1:(strind(1)-2)); % '/media/kuro/DATA1/HCPMEG'
            rMEGPath = rMEGPath(strind(1):end); % '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc'
            anaPath = fullfile([subjectID, '_MEG'],[subjectID, '_MEG_anatomy'],subjectID, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headModelFile = [subjectID '_MEG_anatomy_headmodel.mat'];
            sourceModelFile = [subjectID '_MEG_anatomy_sourcemodel_2d.mat'];
            transformFile = [subjectID '_MEG_anatomy_transform.txt'];
            icaInfo = fullfile(megRootPath,fileparts(rMEGPath),'icaclass',[subjectID '_MEG_3-Restin_icaclass_vs.mat']);
            clear strind;
        elseif isSubjectID
            megRootPath = '/media/kuro/DATA1/HCPMEG'; % Default HCP MEG root dicrectory
            rMEGPath = fullfile([subjectID, '_MEG'],[subjectID, '_MEG_Restin_preproc'],subjectID,'MEG','Restin','rmegpreproc');
            anaPath = fullfile([subjectID, '_MEG'],[subjectID, '_MEG_anatomy'],subjectID, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headModelFile = [subjectID '_MEG_anatomy_headmodel.mat'];
            sourceModelFile = [subjectID '_MEG_anatomy_sourcemodel_2d.mat'];
            transformFile = [subjectID '_MEG_anatomy_transform.txt'];
            rawFile = [subjectID '_MEG_3-Restin_rmegpreproc.mat']; 
            icaInfo = fullfile(megRootPath,fileparts(rMEGPath),'icaclass',[subjectID '_MEG_3-Restin_icaclass_vs.mat']);
        end
    case 'MEG' % This is for future work
        error('pipeline for general MEG data is not ready yet');
    case 'EEG' % This is for future work
        error('pipeline for general EEG data is not ready yet');
    otherwise
        error('unexpected fileType. fileType should be one of ''EEG'', ''MEG'', or ''HCP''.');
end
clear isRawFile isSubjectID subjectID;
%%     PREPROCESSING     %%
%% Load MEG data
clear fileType;
fieldtripData = load(fullfile(megRootPath, rMEGPath, rawFile)); % load HCP rMEG data
fieldtripData = fieldtripData.data;
layoutFile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
try
    D = spm_eeg_load(fullfile(megRootPath,rMEGPath,['fdf_spm_',rawFile]));
catch
    try 
        D = spm_eeg_load(fullfile(megRootPath,rMEGPath,['spm_',rawFile]));
    catch
        D = spm_eeg_ft2spm(fieldtripData, fullfile(megRootPath, rMEGPath,['spm_', rawFile]));
    end

    %% Filters And Downsample %%
    %% High-pass Filter
    S.D                 = D;
    S.band              = 'high';
    S.freq              = 0.1;
    S.prefix            = 'f_';
    D                   = spm_eeg_filter(S);
    clear S;
    %% Downsample
    % S.method - 'resample' (default), 'decimate', 'downsample', 'fft'
    S.D                 = D;
    S.method            = 'fft';
    S.fsample_new       = 200;
    D                   = spm_eeg_downsample(S);
    clear S;
    %% Low-pass Filter
    S.D                 = D;
    S.band              = 'low';
    S.freq              = 100;
    D                   = spm_eeg_filter(S);
    clear S;
    %%
    %% Specify Forward Model
    load(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc
    % read the source and volume conduction model from current dir with
    % outputs of previous pipelines
    fid = fopen(fullfile(megRootPath, anaPath, transformFile), 'rt');
    strFid=fread(fid,[1 inf], '*char');
    eval(strFid);
    fclose(fid);
    fid = [];
    clear fid;
    strFid = [];
    clear strFid;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify sourcemodel, headmodel, sensor and transform if it necessary
    load(fullfile(megRootPath,anaPath,sourceModelFile)); % This file is originally in anatomy directory
    sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
    sourcemodelsubj = sourcemodel2d;

    load(fullfile(megRootPath,anaPath, headModelFile));
    headmodel = ft_convert_units(headmodel, 'cm');

    grad=fieldtripData.grad;

    % These three code lines are for EEG
    % gradBalanced = grad;
    % gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
    % grad=gradBalanced;

    grad = ft_convert_units(grad, 'cm');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare the component data in order for ft_sourceanalysis to be able to
    % swallow it
    channels = comp_class.topolabel;

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
    i = []; headmodel = []; grad = []; ft_default = []; experimentid = []; comp_class = []; channels = []; anaPath = [];
    clear val tlck sourcemodelsubj resultprefix options mixing megpath;
    clear rmeg mainpath i headmodel grad ft_default experimentid comp_class channels anaPath;

    % specify forward model (Bring forward model information 
    D.inv = [];
    D.inv{1}.forward = [];
    D.inv{1}.method='Imaging';
    D.inv{1}.forward.toMNI = transform.bti2spm;
    D.inv{1}.forward.fromMNI = transform.spm2bti;
    D.inv{1}.forward.voltype='Single Shell';
    D.inv{1}.forward.vol = ft_convert_units(job.vol,'m');
    D.inv{1}.forward.modality = 'MEG';
    D.inv{1}.forward.siunits=0;
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

    clear Datareg grid;
    % check meshes and co-registration
    spm_eeg_inv_checkforward(D);
    save(D);
end
%% DCM %%
%% Prepare data ( Set DCM.xY )
DCM = [];
DCM.xY.Dfile = D.fullfile;
DCM.xY.modality = 'MEG';
DCM.name = 'DCM_test';
DCM.val =D.val;
DCM.options.analysis = 'CSD';
DCM.options.model = 'ERP';
DCM.options.spatial = 'ECD';
DCM.options.trials = [1];
DCM.options.Tdcm = [1, 2000];
DCM.options.Fdcm = [1, 100];
DCM.options.onset = 64;
DCM.options.dur = 16;
DCM.options.Nmodes = 8;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
DCM = spm_dcm_erp_data(DCM);
spm_dcm_erp_results(DCM, 'Data');

%% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)
DCM.Lpos = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
DCM.Sname = {'lLP', 'rLP', 'Prec', 'mPFC'};
% DCM.Lpos = [ [-46;19;22] [41;31;30] [-31;21;4] [34;23;1]];
% DCM.Sname = {'lDLPFC','rDLPRF','lVLPRF','rVLPFC'};
Nareas = size(DCM.Lpos,2);
% Spatial Model
DCM = spm_dcm_erp_dipfit_jhs(DCM);

%% Prepare Connectivity Matrices ( Set A, B, C - No need C for 'CSD' )
DCM.xU.name = {}; % If there is some between trial effect, then you should specify xU.
DCM.xU.X = zeros(3,0); % Design Matrix

% Specify connectivity model ( Set as full model )
DCM.A{1} = ones(Nareas, Nareas); % A{1} : Forward
DCM.A{1} = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];
DCM.A{2} = ones(Nareas,Nareas); % A{2} : Backward
DCM.A{2} = DCM.A{1}';
DCM.A{3} = zeros(Nareas,Nareas); % A{3} : Modulatory
DCM.A{3} = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];

DCM.B{1} = zeros(Nareas,Nareas);
%% Invert
DCM = spm_dcm_csd_jhs(DCM);

%% Result Display
Action= {'Spectral data', 'trial-specific effects', 'Input','Transfer functions',...
    'Coupling (A)','Coupling (B)', 'Coupling (C)','Coherence (sources)',...
    'Coherence (channels)','Covariance (sources)','Covariance (channels)',...
    'Cross-spectra (sources)','Cross-spectra (channels)','Dipoles'};
spm_dcm_csd_results(DCM,Action{7});
%% Powermap
S.D                 = D;
S.freqwin           = [8 12];
S.timewin           = [0 2000];
S.freqres           = 1;
S.log               = 0; %1 - yes , 0 - no
Dtf                 = spm_eeg_ft_multitaper_powermap(S);