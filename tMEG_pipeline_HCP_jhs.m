%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP task MEG Source Reconstruction Pipeline Using SPM                   %
%                                                                         %
% HCP task MEG data Simulation & Invert                                   %
% This pipeline code is for HCP dataset                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.01.25 xx:xx - By Junho Son                                     %
%     2018.02.19 16:22 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Path
clear;
% You need to change three variables below for your system.
mPath='/media/kuro/DATA1/HCPMEG';
subjectID='146129';
workingDir = '/home/kuro/Codes/M_code/MEG/HCPtMEG_Simulation_and_Inversion/Test'
%taskType = 'Wrkmem'; % 'Wrkmem',' 


subjectDir = [subjectID,'_MEG'];
tMEGPath = fullfile(subjectDir,[subjectID '_MEG_Wrkmem_preproc'],subjectID,'MEG','Wrkmem','tmegpreproc');
tMEGRawData=[subjectID '_MEG_7-Wrkmem_tmegpreproc_TIM'];
icaInfo = fullfile(mPath,subjectDir,[subjectID '_MEG_Wrkmem_preproc'],subjectID,'MEG','Wrkmem','icaclass',[subjectID,'_MEG_7-Wrkmem_icaclass_vs.mat']);
anaPath = fullfile(subjectDir, [subjectID '_MEG_anatomy'],subjectID, 'MEG','anatomy');
headmodelfile = fullfile(mPath,anaPath,[subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(mPath,anaPath,[subjectID '_MEG_anatomy_sourcemodel_2d']);          
%% Load fieldtrip data
try
    D = spm_eeg_load(fullfile(workingDir,['mapfdfspm_', tMEGRawData]));
    fieldtripData = load(fullfile(mPath, tMEGPath,tMEGRawData));
catch
    fieldtripData = load(fullfile(mPath, tMEGPath, tMEGRawData));
    try
        D = spm_eeg_load(workingDir, ['spm_',tMEGRawData]);
    catch
        D = my_eeg_ft2spm_jhs(fieldtripData,fullfile(mPath, tMEGPath,tMEGRawData), workingDir);
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
    S.timewin = [-500, 1500];
    D = spm_eeg_crop(S);
    S = [];
    S.D = D;
    S.mode = 'reject';
    S.methods.channels = {'all'};
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
end
%% Check ERP signals
labels = D.chanlabels;
for i = 1:9
    figure;
    for j = 1:30
        idx = (i-1)*30+j;
        if idx>size(D,1)
           break; 
        end
        subplot(15,2,j);
        plot(D.time,D(idx,:,1));
        hold on
        plot(D.time,D(idx,:,2));
        plot(D.time,D(idx,:,3));
        hold off
        title(labels(idx));
    end
end
for i = 1:size(D,1)
    r(i)= snr(D(i,:,1), D(i,:,3));
end
ind =find(r>3);
for j = 1:5
    figure;
    for i = 1:10
        idx = (j-1)*10+i;
        if idx>length(ind)
            break
        end
        subplot(5,2,i);
        hold on;
        plot(D.time, D(ind(idx),:,1));
        plot(D.time, D(ind(idx),:,2));
        plot(D.time, D(ind(idx),:,3));
        title(labels(ind(idx)));
        hold off;
    end
end
%% Specify Forward Model
%
%
%
%%%opengl hardware;


% ensure that the time and date of execution are not stored in the provenance information
%%%global ft_default
%%%ft_default.trackcallinfo = 'no';

hcp_read_matlab(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc

% read the source and volume conduction model from current dir with
% outputs of previous pipelines
fid = fopen(fullfile(mPath, anaPath,[subjectID '_MEG_anatomy_transform.txt']), 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify sourcemodel, headmodel, sensor and transform if it necessary
hcp_read_matlab(sourcemodelfile); % This file is originally in anatomy directory
sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
sourcemodelsubj = sourcemodel2d;

hcp_read_matlab(headmodelfile);
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

cfg = [];
cfg.vol = headmodel;
cfg.grid = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;

clear val tlck subjectID sourcemodelsubj resultprefix options mixing megpath;
clear rmeg mainpath i headmodel grad ft_default experimentid comp_class channels anatpath;

% specify forward model (Bring forward model information 
D.inv = [];
D.inv{1}.forward = [];
D.inv{1}.method='Imaging';
D.inv{1}.forward.toMNI = transform.bti2spm;
%
D.inv{1}.forward.fromMNI = transform.spm2bti;
D.inv{1}.forward.voltype='Single Shell';
D.inv{1}.forward.vol = ft_convert_units(cfg.vol,'m');
D.inv{1}.forward.modality = 'MEG';
D.inv{1}.forward.siunits=1;
D.inv{1}.forward.mesh=[];
D.inv{1}.forward.mesh_correction=[];
D.inv{1}.forward.mesh.face = cfg.grid.tri;
grid = ft_convert_units(cfg.grid,'m');
D.inv{1}.forward.mesh.vert = grid.pos;
D.inv{1}.forward.sensors = [];
D.inv{1}.forward.sensors = ft_convert_units(cfg.grad,'m');
D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face;
grid = ft_convert_units(cfg.grid,'mm');
D.inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(D.inv{1}.forward.toMNI,grid.pos); % tess_mni should have mni coordinate
Datareg(1).modality = 'MEG';
D.inv{1}.datareg(1) = Datareg(1);
D.inv{1}.forward.toMNI = [D.inv{1}.forward.toMNI(:, 1:3)*1000 D.inv{1}.forward.toMNI(:,4)];
D.inv{1}.forward.fromMNI = eye(4)/D.inv{1}.forward.toMNI;
clear Datareg grid;
% check meshes and co-registration
spm_eeg_inv_checkforward(D); 
 

%% Invert SPM Example
%  spm_cfg_eeg_inv_invert.m -> function run_inversion(job)
%  Specify job
%

job = [];
job.D = {D}; % job.D is cell type
job.val = 1; %inversion index
job.whatconditions.all = 1;
% Specify Inverting Methods.
% Use 'LOR', 'GS', 'IID', 'EEB' for Methods 'COH', 'MSP(GS)', 'IID', 'EBB'
job.isstandard.custom.invtype='GS'; 
job.isstandard.custom.woi = [0, 600];
job.isstandard.custom.foi = [0, 40]; % simulated data has 10, 20 Hz sources, so 0 to 40 Hz is enough to cover source activity
job.isstandard.custom.hanning = 1;
job.isstandard.custom.priors.priorsmask{1} = '';
job.isstandard.custom.priors.space = 1;
job.isstandard.custom.restrict.locs = [];
job.isstandard.custom.restrict.radius = 32;
job.isstandard.custom.restrict.mask{1}='';
job.modality{1} = 'All';


%
% Run Invert
% spm_cfg_eeg_inv_invert code
%

inverse = [];
if isfield(job.whatconditions, 'condlabel')
    inverse.trials = job.whatconditions.condlabel;
end

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

D = {};

for i = 1:numel(job.D)
    D{i} = spm_eeg_load(job.D{i});
    
    D{i}.val = job.val;
    
    D{i}.con = 1;
    
    if ~isfield(D{i}, 'inv')
        error('Forward model is missing for subject %d.', i);
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning('Duplicating the last forward model for subject %d.', i);
        else
            error('Forward model is missing for subject %d.', i);
        end
    end
    
    D{i}.inv{D{i}.val}.inverse = inverse;
end

D = spm_eeg_invert(D);

if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end

% check invert 
% spm_eeg_invert_display(D{1});
% if you want to see specific condition(eg. condition number 2), use below instead
% D{1}.con = 2;
% spm_eeg_invert_display(D{1});
% check result
% spm_eeg_inv_results(D);