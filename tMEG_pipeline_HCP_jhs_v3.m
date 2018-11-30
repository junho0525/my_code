function DCM = tMEG_pipeline_HCP_jhs_v3(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP Task MEG Source Reconstruction Pipeline Using SPM                   %
%                                                                         %
% This pipeline code is for task HCP dataset                              %
%                                                                         %
% D = tMEG_pipeline_HCP_jhs_v3(cfg)                                       %
%    Output      D     - spm meeg object                                  %
%    Input                                                                %
%        cfg.subjectID     - HCP dataset subject id                       %
%        cfg.tasktype      - 'Wrkmem', 'Motort', 'StoryM'                 %
%        cfg.subjectPath      - HCP Single Subject Data Path              %
%        cfg.savePath      - Save path for output and intermediate files  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.02.20 15:47 - By Junho Son                                     %
%     2018.02.26 19:00 - By Junho Son                                     %
%                Add code for HCP Motor task.                             %
%     2018.02.28 15:16 - By Junho Son                                     %
%                Redefine inputs and options as cfg struct                %
%     2018.03.14 06:44 - By Junho Son                                     %
%                Try not to make segmentation fault                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(cfg, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if ~isfield(cfg, 'taskType'),error('Invalid Task Type!! There is no taskType');end
if (~isfield(cfg, 'subjectPath')||~isdir(cfg.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(cfg, 'savePath')||~isdir(cfg.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(cfg, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end

% These fields is optional. If these fields are not specified, use default.
if ~isfield(cfg, 'sourceLocalization'), cfg.sourceLocalization = 1; end % defulat : 1 
if ~isfield(cfg, 'showERFAll'), cfg.showERFAll = 0; end % defualt : 0
if ~isfield(cfg, 'showForwardResult'), cfg.showForwardResult = 0; end % defualt : 0
if ~isfield(cfg, 'showERFHighSNR'), cfg.showERFHighSNR = 1; end
if ~isfield(cfg, 'showERFSNRThreshold'), cfg.showERFSNRThreshold = 3; end

%% Set Path
switch cfg.taskType
    case 'Wrkmem'
        tMEGPath = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_Wrkmem_preproc'],cfg.subjectID,'MEG','Wrkmem','tmegpreproc');
        tMEGRawData=[cfg.subjectID '_MEG_' num2str(cfg.sessionID) '-Wrkmem_tmegpreproc_TIM'];
        icaInfo = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_Wrkmem_preproc'],cfg.subjectID,'MEG','Wrkmem','icaclass',[cfg.subjectID,'_MEG_' num2str(cfg.sessionID) '-Wrkmem_icaclass_vs.mat']);
    case 'Motort'
        tMEGPath = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_Motort_preproc'],cfg.subjectID,'MEG','Motort','tmegpreproc');
        tMEGRawData=[cfg.subjectID '_MEG_' num2str(cfg.sessionID) '-Motort_tmegpreproc_TFLA'];
        icaInfo = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_Motort_preproc'],cfg.subjectID,'MEG','Motort','icaclass',[cfg.subjectID,'_MEG_' num2str(cfg.sessionID) '-Motort_icaclass_vs.mat']);
    case 'StoryM'
        tMEGPath = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_StoryM_preproc'],cfg.subjectID,'MEG','StoryM','tmegpreproc');
        tMEGRawData=[cfg.subjectID '_MEG_' num2str(cfg.sessionID) '-StoryM_tmegpreproc_TIM'];
        icaInfo = fullfile(cfg.subjectPath,[cfg.subjectID '_MEG_StoryM_preproc'],cfg.subjectID,'MEG','StoryM','icaclass',[cfg.subjectID,'_MEG_' num2str(cfg.sessionID) '-StoryM_icaclass_vs.mat']);
    otherwise
        Error();
end
anaPath = fullfile(cfg.subjectPath, [cfg.subjectID '_MEG_anatomy'],cfg.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[cfg.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[cfg.subjectID '_MEG_anatomy_sourcemodel_2d']);

%% Load fieldtrip data
try
    if cfg.DCM
        switch cfg.DCManalysis
            case 'ERP'
                D = spm_eeg_load(fullfile(cfg.savePath,['mapfdfspm_', tMEGRawData]));
            case 'CSD'
                D = spm_eeg_load(fullfile(cfg.savePath,['apfdfspm_', tMEGRawData]));
            otherwise
        end
    else
        D = spm_eeg_load(fullfile(cfg.savePath,['mapfdfspm_', tMEGRawData]));
    end
    fieldtripData = load(fullfile(tMEGPath,tMEGRawData));
catch
    fieldtripData = load(fullfile(tMEGPath, tMEGRawData));
    try
        D = spm_eeg_load(cfg.savePath, ['spm_',tMEGRawData]);
    catch
        D = my_eeg_ft2spm_jhs(fullfile(tMEGPath,tMEGRawData), cfg.savePath);
    end
    %% 1-20Hz bandpass filter & downsampling
    
    S.D = D;
    S.band = 'high';
    S.freq = 0.9;
    D = spm_eeg_filter(S);
    S = [];
    S.D = D;
    S.method = 'fft';
    S.fsample_new = 200;
    D = spm_eeg_downsample(S);
    S = [];
    S.D = D;
    S.band = 'low';
    S.freq = 20;
    D = spm_eeg_filter(S);
    %% Croping with time window [-500 ms, 1500 ms] , Artefact detaction & Averaging to get ERF
    S = [];
    S.D = D;
    switch cfg.taskType
        case 'Wrkmem'
            S.timewin = [-500, 1500];
        case 'Motort'
            S.timewin = [-300, 1200];
        case 'StoryM'
        otherwise
            error();
    end
    D = spm_eeg_crop(S);
    S = [];
    S.D = D;
    S.mode = 'reject';
    S.methods.channels = {'MEG'};
    S.methods.fun = 'threshchan';
    S.methods.settings.threshold = 3000;%(3000 fT)
    S.methods.settings.excwin = 1000;
    D = spm_eeg_artefact(S);
    S = [];
    S.D = D;
    S.circularise = 0;
    S.robust.ks = 3;
    S.robust.bycondition = 1;
    S.robust.savew = 0;
    S.robust.removebad = 0;
    D = spm_eeg_average(S);
    delete(fullfile(cfg.savePath,['pfdfspm_', tMEGRawData, '.mat']));
    delete(fullfile(cfg.savePath,['pfdfspm_', tMEGRawData, '.dat']));
    delete(fullfile(cfg.savePath,['fdfspm_', tMEGRawData, '.mat']));
    delete(fullfile(cfg.savePath,['fdfspm_', tMEGRawData, '.dat']));
    delete(fullfile(cfg.savePath,['dfspm_', tMEGRawData, '.mat']));
    delete(fullfile(cfg.savePath,['dfspm_', tMEGRawData, '.dat']));
    delete(fullfile(cfg.savePath,['fspm_', tMEGRawData, '.mat']));
    delete(fullfile(cfg.savePath,['fspm_', tMEGRawData, '.dat']));
    delete(fullfile(cfg.savePath,['spm_', tMEGRawData, '.mat']));
    delete(fullfile(cfg.savePath,['spm_', tMEGRawData, '.dat']));
end
%% Check ERP signals
labels = D.chanlabels;
if cfg.saveFig
    mkdir(fullfile(cfg.savePath, 'figure'));
    mkdir(fullfile(cfg.savePath, 'figure', cfg.subjectID));
    mkdir(fullfile(cfg.savePath, 'figure', cfg.subjectID, cfg.taskType));
    savePath = fullfile(cfg.savePath, 'figure', cfg.subjectID, cfg.taskType);
end
if cfg.showERFAll
    for i = 1:9
        figure;
        suptitle(['Deploy ERPs Subject ID : ' cfg.subjectID]);
        for j = 1:30
            idx = (i-1)*30+j;
            if idx>size(D,1)
               break; 
            end
            subplot(15,2,j);
            plot(D.time,D(idx,:,1),'r');
            hold on;
            plot(D.time,D(idx,:,2),'b');
            plot(D.time,D(idx,:,3),'y');
            if strcmp(cfg.taskType ,'Motort')
                plot(D.time, D(idx, :, 4),'g');
                plot(D.time, D(idx, :, 5),'c');
            end
            hold off;
            title(labels(idx));
        end
        if cfg.saveFig
            mkdir(fullfile(savePath, 'allChan'));
            saveas(gcf, fullfile(savePath, 'allChan', ['session_' num2str(cfg.sessionID) '_' num2str(i) '.png']))
        end
        drawnow;
    end
    figure;
    suptitle('All ERFs');
    for i = 1:size(D,3)
        subplot(size(D,3),1,i);
        plot(D.time,D(1,:,i));
        hold on;
        for j = 2:size(D,1)
            plot(D.time,D(j,:,i));
        end
        hold off;
        title(['All ERP - condition_' D.condlist{i}]);
    end
end
if cfg.showERFHighSNR
    for i = 1:size(D,1)
        switch cfg.taskType
            case 'Wrkmem'
                r(i)= snr(D(i,:,1), D(i,:,3));
            case 'Motort'
                r(i)= snr(D(i,:,1), D(i,:,5));
            otherwise
        end
    end
    ind =find(r>cfg.showERFSNRThreshold);
    for j = 1:5
        idx = (j-1)*10;
        if idx>=length(ind)||idx>=cfg.showERFSNRMaxNum
            break;
        end
        figure;
        for i = 1:10
            idx = (j-1)*10+i;
            if idx>length(ind)||idx>cfg.showERFSNRMaxNum
                break;
            end
            subplot(5,2,i);
            hold on;
            plot(D.time, D(ind(idx),:,1),'r');
            plot(D.time, D(ind(idx),:,2),'b');
            if strcmp(cfg.taskType ,'Motort')
                plot(D.time, D(ind(idx),:,3),'g');
                plot(D.time, D(idx, :, 4),'c');
                plot(D.time, D(idx, :, 5),'black');
            elseif strcmp(cfg.taskType ,'Wrkmem')
                plot(D.time, D(ind(idx),:,3),'black');
            end
            title(labels(ind(idx)));
            hold off;
        end
        suptitle(['Channels that show high SNR ([Condition1 Signal]/[Fixation Signal]) subject ID = ' cfg.subjectID]);
        if cfg.saveFig
            mkdir(fullfile(savePath, 'highSNR'));
            saveas(gcf ,fullfile(savePath, 'highSNR',['session_' num2str(cfg.sessionID) '_' num2str(j) '.png']));
        end
        drawnow;
        close(gcf);
    end
    labels = [];
    clear labels;
end
%% Source Localization
if cfg.sourceLocalization
    %% Specify Forward Model
    load(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc
    % read the source and volume conduction model from current dir with
    % outputs of previous pipelines
    fid = fopen(fullfile(anaPath,[cfg.subjectID '_MEG_anatomy_transform.txt']), 'rt');
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
    save(fullfile(cfg.savePath,D.fname),D);
    
    clear Datareg grid;
    % check meshes and co-registration
    if cfg.showForwardResult; spm_eeg_inv_checkforward(D); end
    %% Invert SPM Example
    %  spm_cfg_eeg_inv_invert.m -> function run_inversion(job)
    
    % Specify job
    job = [];
    job.D = {D}; % job.D is cell type
    job.val = 1; %inversion index
    job.whatconditions.all = 1;
    % Specify Inverting Methods.
    % Use 'LOR', 'GS', 'IID', 'EEB' for Methods 'COH', 'MSP(GS)', 'IID', 'EBB'
    job.isstandard.custom.invtype='GS'; 
    job.isstandard.custom.woi = [-200, 1000];
    job.isstandard.custom.foi = [0, 40]; % simulated data has 10, 20 Hz sources, so 0 to 40 Hz is enough to cover source activity
    job.isstandard.custom.hanning = 1;
    job.isstandard.custom.priors.priorsmask{1} = '';
    job.isstandard.custom.priors.space = 1;
    job.isstandard.custom.restrict.locs = [];
    job.isstandard.custom.restrict.radius = 32;
    job.isstandard.custom.restrict.mask{1}='';
    job.modality{1} = 'MEG';

    % Run Invert (spm_cfg_eeg_inv_invert code)
    inverse = [];
    if isfield(job.whatconditions, 'condlabel'),inverse.trials = job.whatconditions.condlabel;end
    if isfield(job.isstandard, 'custom')
        inverse.type = job.isstandard.custom.invtype;
        inverse.woi  = fix([max(min(job.isstandard.custom.woi), 1000*D.time(1)) min(max(job.isstandard.custom.woi), 1000*D.time(end))]);
        inverse.Han  = job.isstandard.custom.hanning;
        inverse.lpf  =  fix(min(job.isstandard.custom.foi));
        inverse.hpf  =  fix(max(job.isstandard.custom.foi));
        
        P = char(job.isstandard.custom.priors.priorsmask);
        if ~isempty(P)        
            [p,f,e] = fileparts(P);
            switch lower(e)
                case '.gii'
                    g = gifti(P);
                    inverse.pQ = cell(1,size(g.cdata,2));
                    for i=1:size(g.cdata,2)
                        inverse.pQ{i} = double(g.cdata(:,i));
                    end
                case '.mat'
                    load(P);
                    inverse.pQ = pQ;
                case {'.img', '.nii'}
                    S.D = D;
                    S.fmri = P;
                    S.space = job.isstandard.custom.priors.space;
                    D = spm_eeg_inv_fmripriors(S);
                    inverse.fmri = D.inv{D.val}.inverse.fmri;
                    load(inverse.fmri.priors);
                    inverse.pQ = pQ;
                otherwise
                    error('Unknown file type.');
            end
        end

        if ~isempty(job.isstandard.custom.restrict.locs)
            inverse.xyz = job.isstandard.custom.restrict.locs;
            inverse.rad = job.isstandard.custom.restrict.radius;
        end

        P = char(job.isstandard.custom.restrict.mask);
        if ~isempty(P)
            inverse.mask = P;
        end
    end

    [mod, list] = modality(D, 1, 1);
    if strcmp(job.modality{1}, 'All')
        inverse.modality  = list;
    else
        inverse.modality  = intersect(list, job.modality);
    end

    if numel(inverse.modality) == 1
        inverse.modality = inverse.modality{1};
    end
    D.val = job.val;
    D.con = 1;
    if ~isfield(D, 'inv'), error('Forward model is missing for subject %d.', cfg.subjectID);
    elseif numel(D.inv)<D.val || ~isfield(D.inv{D.val}, 'forward'), 
        if D.val>1 && isfield(D.inv{D.val-1}, 'forward')
            D.inv{D.val} = D.inv{D.val-1};
            warning('Duplicating the last forward model for subject %d.', cfg.subjectID);
        else
            error('Forward model is missing for subject %d.', cfg.subjectID);
        end
    end
    D.inv{D.val}.inverse = inverse;
    D = spm_eeg_invert(D);
    drawnow;
    if iscell(D)
        D = D{1};
    end
    save(fullfile(cfg.savePath,D.fname),D);
    if cfg.saveFig
        mkdir(fullfile(savePath, 'SourceLoc'));
        switch cfg.taskType
            case 'Wrkmem'
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_0back.png']));
                D.con=2;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_2back.png']));
                D.con=3;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_fix.png']));
            case 'Motort'
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_LH.png']));
                f= gcf;
                pause(0.1);
                close(f);
                waitfor(f);
                D.con=2;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_LF.png']));
                f= gcf;
                pause(0.1);
                close(f);
                waitfor(f);
                D.con=3;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_RH.png']));
                f= gcf;
                pause(0.1);
                close(f);
                waitfor(f);
                D.con=4;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_RF.png']));
                f= gcf;
                pause(0.1);
                close(f);
                waitfor(f);
                D.con=5;
                spm_eeg_invert_display(D);
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source_fix.png']));
                f= gcf;
                pause(0.1);
                close(f);
                waitfor(f);
            otherwise
                saveas(gcf ,fullfile(savePath, 'SourceLoc', [cfg.subjectID '_' num2str(cfg.sessionID) '_Source.png']));
        end
    end
    if ~cfg.DCM
        cfg = []; fieldtripData = []; headmodelfile = []; icaInfo =[]; idx = []; ind = []; inverse = []; j = []; job = [];
        list = []; mod = []; P = []; r = []; resultprefix = []; savePath = []; sourcemodel2d = []; sourcemodelfile = [];
        tMEGPath = []; tMEGRawData = []; transform=[];
        clear cfg fieldtirpData headmodelfile icaInfo idx ind inverse j job list mod P r resultprefix savePath sourcemodel2d sourcemodelfile tMEGPath;
        clear tMEGRawData transform;
        cfg.DCM = 0;
        %% D = [];
    end
end
%% DCM
if cfg.DCM
    %% Check forward
    try
        clear D;
        switch cfg.DCManalysis
            case 'ERP'
                D = spm_eeg_load(fullfile(cfg.savePath,['mapfdfspm_', tMEGRawData]));
            case 'CSD'
                D = spm_eeg_load(fullfile(cfg.savePath,['mapfdfspm_', tMEGRawData]));
        end
        if isfield(D,'inv') &&isfield(D.inv{1},'forward')&&isfield(D.inv{1},'datareg')&&isfield(D.inv{1},'mesh')
        else
            %% Specify Forward Model
            load(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc
            % read the source and volume conduction model from current dir with
            % outputs of previous pipelines
            fid = fopen(fullfile(anaPath,[cfg.subjectID '_MEG_anatomy_transform.txt']), 'rt');
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
            save(fullfile(cfg.savePath,D.fname),'D');
            clear Datareg grid;
            % check meshes and co-registration
            if cfg.showForwardResult; spm_eeg_inv_checkforward(D); end
            switch cfg.DCManalysis
                case 'ERP'
                    save(D.fullfile,'D');
                case 'CSD'
                    save(D.fullfile,'D');
            end
        end
    catch
        error();
    end   
    %% Specify DCM
    DCM.name = fullfile(cfg.savePath,strcat('DCM_',tMEGRawData));
    DCM.xY.Dfile          = D.fullfile;
    
    dt = 5; % ms; sampling rate = 200 Hz
    DCM.options.analysis  = cfg.DCManalysis;
    DCM.options.spatial   = 'ECD';
    DCM.options.model     = 'ERP';
    try
        DCM.options.Tdcm = cfg.DCMTimeWin;
    catch
        switch cfg.taskType
            case 'Wrkmem'
                DCM.options.Tdcm   = [0, 500];  % [start end] time window in ms
                % DCM.options.onset  = 0;     % onset time of the first stimulus (in ms)
                % DCM.options.dur    = 150;      % duration
            case 'Motort'
                DCM.options.Tdcm = [-300, 1200];  % [start end] time window in ms
    %             DCM.options.onset  = 0;     % onset time of the first stimulus (in ms)
    %             DCM.options.dur    = 150;      % duration
            case 'StoryM'
            otherwise
                error();
        end
    end
    try 
        DCM.options.trials = cfg.DCMTrial;
    catch
        switch cfg.taskType
            case 'Wrkmem'
                DCM.options.trials = [1 2 3];
            case 'Mortort'
            otherwise
        end
    end
    try
        DCM.options.Fdcm = cfg.DCMFreqWin;
    catch
        DCM.options.Fdcm      = [1,60];     % [start end] Frequency window in Hz
    end
    DCM.options.Nmodes    = 8;      % # of feature
    DCM.options.DATA      = 1;      % y/n preprocessed data
    DCM.options.D         = 1;      % Down-sampling interval
    DCM.options.onset     = 0;
    DCM.options.dur       = 150;
    DCM.options.h         = 1;
    DCM.options.lock      = 1;
    DCM.options.multiC    = 0;
    DCM.options.symmetry  = 0;
    DCM.options.Rft       = 5;
    
    
    DCM = spm_dcm_erp_data(DCM);
    spm_dcm_erp_results(DCM, 'Data');
    %% Location Priors for dipoles
    try 
        DCM.Lpos  = cfg.DCMROI;
        Nareas    = size(DCM.Lpos,2);
        try 
            DCM.Sname = cfg.DCMROIname; 
        catch
            DCM.Sname = cell(1,Nareas);
            for i =1:Nareas
                DCM.Sname{i} = ['roi_',i];
            end
        end
    catch
        DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
        DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};
        Nareas    = size(DCM.Lpos,2);
    end
    
    % Spatial Model
    DCM = spm_dcm_erp_dipfit_jhs(DCM);
    
    %% Prepare Connectivity Matrices ( A, B, C )
    if isfield(cfg, 'DCMAMatrix')
        DCM.A = cfg.DCMAMatrix;
    else
        % By Default, make full model
        DCM.A{1} = tril(ones(Nareas,Nareas))-eye(Nareas); % A{1} : forward
        DCM.A{2} = triu(ones(Nareas,Nareas))-eye; % A{2} : backward
        DCM.A{3} = ones(Nareas,Nareas); % A{3} : lateral
    end
    if isfield(cfg, 'DCMBMatrix')
        DCM.B = cfg.DCMBMatrix;
    else
        if isfield(cfg,'DCMeffect')
            for i = 1:size(cfg.DCMeffect,2)
                DCM.B{i} = ones(Nareas,Nareas);
            end
        else
            DCM.B{1} = zeros(Nareas,Nareas);
        end
    end
    if isfield(cfg,'DCMCMatrix')
        DCM.C = cfg.DCMCMatrix;
    else
        DCM.C = [0;0;0;0;1];
    end
    if isfield(cfg,'DCMeffect')
        DCM.xU.X = cfg.DCMeffect;
    else
        DCM.xU.X = zeros(length(D.condlist),0);
    end
    try
        DCM.xU.name = cfg.DCMeffectName;
    catch
    end
    %% Invert
    switch cfg.DCManalysis
        case 'CSD'
            DCM = spm_dcm_csd_jhs(DCM);
        case 'ERP'
            DCM = spm_dcm_erp_jhs(DCM);
        otherwise
    end
end