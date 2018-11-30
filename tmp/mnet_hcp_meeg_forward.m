function [D,fieldtripData] = mnet_hcp_meeg_forward(config,D,fieldtripData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare MEEG Forward Model for HCP dataset by using SPM                 %
%                                                                         %
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
icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-Restin_icaclass_vs.mat']);
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
%% Load MEEG Datas
if nargin==1 || isempty(D)||isempty(fieldtripData)
    try
        D = spm_eeg_load(fullfile(config.savePath,['affdspm_', rMEGRawData]));
        fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
    catch
        error(['There is no preprocessed file!! filename : ' fullfile(config.savePath,['affdspm_', rMEGRawData])]);
    end
end
%% Check forward
if isfield(D,'inv') &&isfield(D.inv{1},'forward')&&isfield(D.inv{1},'datareg')&&isfield(D.inv{1},'mesh')
else
    %% Specify Forward Model
    load(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc
    % read the source and volume conduction model from current dir with
    % outputs of previous pipelines
    fid = fopen(fullfile(anaPath,[config.subjectID '_MEG_anatomy_transform.txt']), 'rt');
    strFid=fread(fid,[1 inf], '*char');
    eval(strFid);
    fclose(fid);
    fid = [];
    clear fid;
    strFid = [];
    clear strFid;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify sourcemodel, headmodel, sensor and transform if it necessary
    load(sourcemodelfile); % This file is originally in anatomy directory
    sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
    sourcemodelsubj = sourcemodel2d;

    load(headmodelfile);
    headmodel = ft_convert_units(headmodel, 'cm');

    grad=fieldtripData.data.grad;

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
    save(fullfile(config.savePath,D.fname),'D');
    clear Datareg grid;
    % check meshes and co-registration
    % if config.showForwardResult; spm_eeg_inv_checkforward(D); end
    save(fullfile(config.savePath,D.fname),'D');
end
