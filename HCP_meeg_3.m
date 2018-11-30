%% HCP MEG Data를 DCM으로 돌리기 위한 routine

%% HCP MEG Data는 이미 어느정도 Preprocessing이 되어있다.
%  HCP MEG Processing Pipeline을 통해, datacheck, baddata, icaclass를 통한 data
%  reduction이 이루어져 있는 상태이다.
%  Resting state analysis의 경우, rmegpreproc을, 사용하자.


%%  Load HCP MEG DATA 
%   나중에 배치로 돌리기위해 수정하자.
close all;
clear;
cd ('/home/kuro/Codes/M_code/Simul/HCP');
% datasetPath = '/media/kuro/DATA1/HCPMEG';
% workingPath = '/home/kuro/Codes/M_code/HCP_meeg';
% %datasetPath = '/home/kuro/Codes/M_code/SPM_Examples/EEG_Single/';
% cd (datasetPath);
% apath = '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc';
filename = '100307_MEG_3-Restin_rmegpreproc.mat';
%filename = 'subject1.bdf';


fieldtripData = load(filename); % MEG Data는 field trip data로 저장되어있다.
%fieldtripData = spm_eeg_convert([datasetPath,filename]);

if isfield( fieldtripData,'time')
    data = spm_eeg_ft2spm(fieldtripData, 'spm_MEG');
else
    data = spm_eeg_ft2spm(fieldtripData.data, 'spm_MEG');
end
%   [f,p]        = uigetfile('*.mat','please select data file');
%   file이 안열릴때는 위의 코드를 이용해 유저가 직접 파일을 열수있도록 하자.
%%  Resting State Data Vectorize
%   왠지는 모르겠으나, HCP데이터는 하나의 time-series data를 여러 trial로 쪼갠다(resting임에도.)
%   따라서, 이를 하나의 time-series로서 재구성하는 과정.


%%  SPM MEG Preprocessing
%   SPM EEG 예제의 Preprocessing 절차를 그대로 따라간다.
%   High-pass filter -> Downsampling -> Low-pass filter
%    (Cutoff = 0.1)   (down to 200Hz)    (Cutoff = 30)
%   다만 Slow wave의 경우, 40Hz이상을 보기도 하기때문에, 이경우 Low-pass filter의 
%   Cutoff를 더 높게 조정해주면 된다.
%   filter의 경우, spm_eeg_filter를 사용한다.
%       FORMAT D = spm_eeg_filter(S)
%
%       S           - input structure
%        Fields of S:
%         S.D       - MEEG object or filename of M/EEG mat-file
%
%         S.band    - filterband [low|high|bandpass|stop]
%         S.freq    - cutoff frequency(-ies) [Hz]
%
%        Optional fields:
%         S.type    - filter type [default: 'butterworth']
%                         'butterworth': Butterworth IIR filter
%                         'fir':         FIR filter (using MATLAB fir1 function)
%         S.order   - filter order [default: 5 for Butterworth]
%         S.dir     - filter direction [default: 'twopass']
%                         'onepass'         forward filter only
%                         'onepass-reverse' reverse filter only, i.e. backward in time
%                         'twopass'         zero-phase forward and reverse filter
%         S.prefix  - prefix for the output file [default: 'f']
%
%        D           - MEEG object (also written to disk)

%%
%   High-pass Filter
S.D = data;
S.band = 'high';
S.freq = 0.1;
S.prefix = 'f_';
data = spm_eeg_filter(S);
clear S;

%%
%   Downsampling
S.D = data;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
data = spm_eeg_downsample(S);
clear S;

%%
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = data;
S.band = 'low';
S.freq = 30;
data = spm_eeg_filter(S);
clear S;

%% Specify cfg and some configurations
%mnet_hcp_meeg_defaults;
opengl hardware;
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
return_path = pwd;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

%experimentid = 'MEG';
%scanid = 'Resing';
subjectid = '100307';
experimentid=[subjectid '_MEG'];
mainpath=['/media/kuro/DATA1/HCPMEG/' experimentid];
megpath=fullfile(mainpath,[experimentid '_Restin_preproc'],subjectid,'MEG','Restin');
anatpath=fullfile(mainpath,[experimentid '_anatomy'],subjectid,'MEG','anatomy');
resultprefix=fullfile(megpath,'icaclass',[subjectid '_MEG_3-Restin']);


cd(megpath)
hcp_read_matlab([resultprefix '_icaclass_vs.mat'])

% read the source and volume conduction model from current dir with
% outputs of previous pipelines



cd(anatpath);
fid = fopen('100307_MEG_anatomy_transform.txt', 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify sourcemodel, headmodel, sensor and transform if it necessary
hcp_read_matlab([subjectid '_MEG_anatomy_sourcemodel_2d.mat']);
sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
%sourcemodel2d = ft_transform_geometry(transform.bti2spm, sourcemodel2d);
sourcemodelsubj = sourcemodel2d;

%sourcemodel2d.inside = 1:size(sourcemodel2d.pos,1);
%sourcemodel2d.outside = [];

hcp_read_matlab(sprintf('%s.mat', [subjectid '_MEG_anatomy_headmodel']));
headmodel = ft_convert_units(headmodel, 'mm');
%headmodel = ft_transform_geometry(transform.bti2spm, headmodel);
headmodel = ft_convert_units(headmodel, 'cm');

cd(megpath)

grad=fieldtripData.data.grad;
% gradBalanced = grad;%
% gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');%
% grad=gradBalanced;%
grad = ft_convert_units(grad, 'mm');
%grad = ft_transform_geometry(transform.bti2spm, grad);
grad = ft_convert_units(grad, 'cm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the component data in order for ft_sourceanalysis to be able to
% swallow it
mixing   = comp_class.topo;
channels = comp_class.topolabel;
% normalisation of the topographies
for i = 1:size(mixing, 2)
  val(i) = 0.01*max(abs(mixing(:, i)));
  mixing(:, i) = mixing(:, i)/val(i);
end

% create a 'timelock' structure
tlck = [];
tlck.label = channels;
tlck.cov = eye(numel(tlck.label)); 
tlck.time=1;
tlck.grad = grad;
tlck.dimord = 'chan_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the forward solution

cfg = [];
cfg.vol = headmodel;
cfg.grid = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;

clear val tlck subjectid sourcemodelsubj resultprefix options mixing megpath;
clear rmeg mainpath i headmodel grad ft_default experimentid comp_class channels anatpath;

%% specify forward model (이미 HCP에서 계산해놓은 Forward Model을 가져오는 과정)
cd (return_path);
D = data;
clear data;
cd('..');
fid = fopen('100307_MEG_anatomy_transform.txt', 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;
D.inv = [];
D.inv{1}.forward = [];
D.inv{1}.method='Imaging';
D.inv{1}.forward.toMNI = transform.bti2spm;
D.inv{1}.forward.fromMNI = transform.spm2bti;
% D.inv{1}.forward.toMNI = eye(4,4);
% D.inv{1}.forward.fromMNI = eye(4,4);

D.inv{1}.forward.voltype='Single Shell';
D.inv{1}.forward.vol = ft_convert_units(cfg.vol,'m'); % 

D.inv{1}.forward.modality = 'MEG';
D.inv{1}.forward.siunits=1;
D.inv{1}.forward.mesh=[];
D.inv{1}.forward.mesh_correction=[];


cfg.grid = ft_convert_units(cfg.grid, 'mm');
D.inv{1}.forward.mesh.face = cfg.grid.tri;
D.inv{1}.forward.mesh.vert = cfg.grid.pos;


D.inv{1}.forward.sensors = [];
D.inv{1}.forward.sensors = ft_convert_units(cfg.grad,'m'); %



D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face; %% mesh.tess_mni는 mm로 되어있어야 하며
%D.inv{1}.mesh.tess_mni.vert = D.inv{1}.forward.mesh.vert;
D.inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(D.inv{1}.forward.toMNI, D.inv{1}.forward.mesh.vert);% mm

%D.inv{1}.mesh.tess_mni.face = sourcemodel2d.tri; %% mesh.tess_mni는 mm로 되어있어야 하며
%D.inv{1}.mesh.tess_mni.vert = sourcemodel2d.pos;

%D.inv{1}.mesh.tess_mni.vert = [D.inv{1}.forward.mesh.vert(:,2),D.inv{1}.forward.mesh.vert(:,1),D.inv{1}.forward.mesh.vert(:,3)];%
D.inv{1}.forward.mesh.vert = D.inv{1}.forward.mesh.vert./1000; %% mesh는 m로 되어 있어야 제대로 기능한다.
Datareg(1).modality = 'MEG';
D.inv{1}.datareg(1) = Datareg(1);



clear Datareg;
% check meshes and co-registration
% spm_eeg_inv_checkforward(D); 

%% Simulation Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% specify job for spm_cfg_eeg_inv_simulate.m, function run_simulation(job)
cd ('HCP');
job = [];
job.D = [];
job.D{1} = 'fdf_spm_MEG';
job.val = 1;
job.prefix = 'sim_';
job.whatconditions.all =1;
job.isinversion.setsources.woi = [100 400];

job.isinversion.setsources.isSin.foi = [10 ; 20];% alpha, beta
%job.isinversion.setsources.isSin.foi = ones( size(D.inv{1}.mesh.tess_mni.vert,1),1)*20;

%job.isinversion.setsources.dipmom=[5,10;10,5;10,3];
job.isinversion.setsources.dipmom=[10 10; 10 5];
%job.isinversion.setsources.dipmom=ones(size(D.inv{1}.mesh.tess_mni.vert,1),2)*10;

%job.isinversion.setsources.locs = [-6 28 34; -6 -52 40; 87 74 86]; %ACC, PCC, Precentral gyrus
%job.isinversion.setsources.locs = [-6 28 34; -6 -52 40; -46 -70 36 ]; %ACC, PCC, LLP
job.isinversion.setsources.locs = [53.7203     -25.7363       9.3949; -52.8484     -25.7363       9.3949 ]; %ACC, PCC, LLP
%job.isinversion.setsources.locs = D.inv{1}.mesh.tess_mni.vert;
%job.isinversion.setsources.locs = [D.inv{1}.mesh.tess_mni.vert(:,2),D.inv{1}.mesh.tess_mni.vert(:,1),D.inv{1}.mesh.tess_mni.vert(:,3)];
job.isSNR.setSNR = 0;

%job.isinversion.setsources.locs = spm_eeg_inv_transform_points(D.inv{1}.forward.fromMNI,job.isinversion.setsources.locs);

%%
% modified code of spm_cfg_eeg_inv_simulate.m
trialind=[];
if isfield(job.whatconditions, 'condlabel')
    trialind =D.indtrial( job.whatconditions.condlabel);
    if isempty(trialind),
        error('No such condition found');
    end;
end
if numel(job.D)>1,
    error('Simulation routine only meant to replace data for single subjects');
end;



if size(modality(D),1)>1,
    error('only suitable for single modality data at the moment');
end;

save(D);
D = {};


for i = 1:numel(job.D) %% only set up for one subject at the moment but leaving this for the future
    D{i} = spm_eeg_load(job.D{i});
    D{i}.val = job.val;
    
    D{i}.con = 1;
    if ~isfield(D{i}, 'inv')
        error(sprintf('Forward model is missing for subject %d', i));
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning(sprintf('Duplicating the last forward model for subject %d', i));
        else
            error(sprintf('Forward model is missing for subject %d', i));
        end
    end
    
end; % for i


if isfield(job.isSNR,'whitenoise'),
    whitenoisefT=job.isSNR.whitenoise; %% internal units as femto Tesla
    SNRdB=[];
else
    SNRdB=job.isSNR.setSNR;
    whitenoisefT=[];
end;

if isfield(job.isinversion,'setsources'), %% defining individual sources (이 블록이 실행되야함)
    
    %%%
    Nsources=size(job.isinversion.setsources.locs,1)
    
    if (size(job.isinversion.setsources.dipmom,1)~=Nsources),
        error('Number of locations must equal number of moments specified');
    end;
    
    mnimesh=[]; %% using mesh defined in forward model at the moment
    
    
    ormni=[]; %% dipoles will get orientations from cortical mesh
    dipmom=job.isinversion.setsources.dipmom;
    if size(job.isinversion.setsources.dipmom,2)==3, %% dipole orientation is specified 3d인경우
        disp('Simulating dipoles without reference to cortical surface');
        for j=1:size(dipmom,1),
            ormni(j,:)=dipmom(j,:)./sqrt(dot(dipmom(j,:),dipmom(j,:))); %% unit orientation in MNI space
            nAmdipmom(j)=sqrt(dot(dipmom(j,:),dipmom(j,:))); % magnitude of dipole
        end;
    else %% only one moment parameter given
        nAmdipmom=dipmom(:,1); %% total momnent in nAm
        dipfwhm=[];
        if size(dipmom,1)==2,
            dipfwhm=dipmom(:,2); %% fhwm in mm
        end;
    end;
    
    woi=job.isinversion.setsources.woi./1000;
    timeind = intersect(find(D{1}.time>woi(1)),find(D{1}.time<=woi(2)));
    simsignal = zeros(Nsources,length(timeind));
    
    if isfield(job.isinversion.setsources.isSin,'fromfile'),
        % Simulate orthogonal Gaussian signals
        filename=cell2mat(job.isinversion.setsources.isSin.fromfile);
        a=load(filename,'s');
        
        if size(a.s,1)~=Nsources,
            error('size of simulated data in file does not match number of sources');
        end;
        
        if size(a.s,2)~=length(timeind),
            warning(sprintf('sim signal from file is not same length as time window (%d vs %d samples) truncating or padding with zeros',size(a.s,2),length(timeind)));
        end;
        usesamples=1:min(length(timeind),size(a.s,2));
        simsignal(:,usesamples)=a.s(:,usesamples);
        
    end; % if isfield fband
    
    if isfield(job.isinversion.setsources.isSin,'fband'),
        % Simulate orthogonal Gaussian signals
        
        simsignal=randn(Nsources,length(timeind));
        % filter to bandwidth
        simsignal=ft_preproc_lowpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(2),2);
        simsignal=ft_preproc_highpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(1),2);
        [u,s,v]=svd(simsignal*simsignal');
        simsignal=u*simsignal; %% orthogonalise all signals
    end; % if isfield fband
    
    %  simsignal=simsignal.*repmat(nAmdipmom,1,size(simsignal,2)); %% now scale by moment
    
    if isfield(job.isinversion.setsources.isSin,'foi'),
        % simulate sinusoids
        sinfreq=job.isinversion.setsources.isSin.foi;
        % Create the waveform for each source
        
        for j=1:Nsources                % For each source
            simsignal(j,:)=sin((D{1}.time(timeind)- D{1}.time(min(timeind)))*sinfreq(j)*2*pi);
        end; % for j
        
        
    end; %% if isfield foi
    
    simsignal=simsignal./repmat(std(simsignal'),size(simsignal,2),1)'; %% Set sim signal to have unit variance
    %[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.isinversion.setsources.locs,simsignal,woi,whitenoisefT,SNRdB,trialind,mnimesh,SmthInit);
    figure();
    plot(D{1}.time(timeind),simsignal);
    xlabel('time');
    ylabel('Normalized amplitude (to be scaled by dip moment later)');
    legend(num2str([1:Nsources]'));
    [D,meshsourceind]=spm_eeg_simulate(D,job.prefix,job.isinversion.setsources.locs,simsignal,ormni,woi,whitenoisefT,SNRdB,trialind,mnimesh,dipfwhm,nAmdipmom);
    
else %% simulate sources based on inversion
    if ~isfield(D{i}.inv{job.val},'inverse'),
        error('There is no solution defined for these data at that index');
    end;
    
    [D]=spm_eeg_simulate_frominv(D,job.prefix,job.val,whitenoisefT,SNRdB,trialind);
end;



if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end


fullname=[D{1}.path filesep D{1}.fname];
out.D = {fullname};


%% Inversion Code
%  spm_cfg_eeg_inv_invert.m -> function run_inversion(job)
%  run_inversion을 위한 job 설정.
%  simulated data를 inversion해 simulation에서 넣어준 Source location을 잘 추정하는지 보자.

job = [];
job.D = D;
job.val = 1; %inversion index
job.whatconditions.all = 1;
% 'COH', 'MSP(GS)', 'IID', 'EBB'를 각각 'LOR', 'GS', 'IID', 'EEB'로 설정해주면 된다.
% 이중 하나를 골라 설정.
job.isstandard.custom.invtype='GS'; 
job.isstandard.custom.woi = [0, 600];
job.isstandard.custom.foi = [0, 80]; % simulated data는 10, 20Hz를 사용하므로 이정도면 된다.
job.isstandard.custom.hanning = 1;
job.isstandard.custom.priors.priorsmask{1} = '';
job.isstandard.custom.priors.space = 1;
job.isstandard.custom.restrict.locs = [];
job.isstandard.custom.restrict.radius = 32;
job.isstandard.custom.restrict.mask{1}='';
job.modality{1} = 'All';


%% spm_cfg_eeg_inv_invert code
D = D{1};
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


