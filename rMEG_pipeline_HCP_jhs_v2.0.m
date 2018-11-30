%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This Is M/EEG Resting State Source Localization Pipeline Code      %
%                           By Using Fieldtrip                            %
%                                                                         %
% We need five inputs                                                     %
% You need to specify these parameters.( All of them are string )         %
%     filetype        - File type               ('EEG', 'MEG' or 'HCP')   %
%     rawfile         - M/EEG raw data filename (eg. Subject2Rest.ds)     %
%     subjectid       - Subject ID for HCP data (eg. '100307')            %
%     sourcemodelfile - Source model filename   (eg. sourcemodel4k.mat)   %
%     headmodelfile   - Head model filename     (eg. headmodel.mat)       %
%     layoutfile      - Channel layout file     (eg. 4D248_helmet.mat)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.02.09 17:39 - By Junho Son                                     %
%     2018.03.30 22:04 - By Junho Son                                     %
%                Add code for FC(coh) and full voxel network analysis     %
%     2018.03.31 01:00 - By Junho Son                                     %
%                Add Seed based FC                                        %
%     2018.04.02 20:10 - By Junho Son                                     %
%                Show FC(coh) on cortical surface                         %
%     2018.05.24 20:31 - By Junho Son                                     %
%                Add source estimation with 'dics' and split trials into  %
%                two groups                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SET DATA PATH    %%
%% Set Library Path
addpath('/home/kuro/Codes/M_code/my_code');
use_my_library('spm',0);
use_my_library('ft',1);
%% Check Input Parameters
%  If you want to use example datapath, please excute this section twice.
if ~exist('filetype')
   warning('You should specify file type.');
   filetype = 'HCP';
end
if exist('rawfile')
    israwfile=1;
else
    israwfile =0;
end
if exist('subjectid')
    issubjectid = 1;
else
    issubjectid = 0;
end

% distract path from rawfile
if ~israwfile&&~issubjectid
    warning('There is no rawfile. This code use default dataset (HCP-subject 100307 resting data)');
    % rawfile = 'D:\HCPDATA\100307_MEG\100307_MEG_Restin_preproc\100307\MEG\Restin\rmegpreproc\100307_MEG_3-Restin_rmegpreproc';
    rawfile = '/media/kuro/DATA1/HCPMEG/100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat';
end
switch filetype
    case 'HCP'
        if israwfile
            [rmegpath rawfile]    = fileparts(rawfile);
            subjectid = strtok(rawfile, '_');
            strind = strfind(rmegpath, subjectid);
            
            megrootpath = rmegpath(1:(strind(1)-2)); % '/media/kuro/DATA1/HCPMEG'
            rmegpath = rmegpath(strind(1):end); % '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc'
            anapath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_anatomy'],subjectid, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headmodelfile = [subjectid '_MEG_anatomy_headmodel.mat'];
            sourcemodelfile = [subjectid '_MEG_anatomy_sourcemodel_2d.mat'];
            coordtransformfile = [subjectid '_MEG_anatomy_transform.txt'];
            
            clear strind;
        elseif issubjectid
            megrootpath = '/media/kuro/DATA1/HCPMEG'; % Default HCP MEG root dicrectory
            rmegpath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_Restin_preproc'],subjectid,'MEG','Restin','rmegpreproc');
            anapath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_anatomy'],subjectid, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headmodelfile = [subjectid '_MEG_anatomy_headmodel.mat'];
            sourcemodelfile = [subjectid '_MEG_anatomy_sourcemodel_2d.mat'];
            coordtransformfile = [subjectid '_MEG_anatomy_transform.txt'];
            rawfile = [subfectid '_MEG_3-Restin_rmegpreproc.mat']; 
        end
    case 'MEG' % This is for future work
        error('pipeline for general MEG data is not ready yet');
    case 'EEG' % This is for future work
        error('pipeline for general EEG data is not ready yet');
    otherwise
        error('unexpected filetype. filetype should be one of ''EEG'', ''MEG'', or ''HCP''.');
end
clear israwfile issubjectid subjectid;
%%     PREPROCESSING     %%
%% Load MEG data
if strcmp('HCP', upper(filetype))
    clear filetype;
    load(fullfile(megrootpath, rmegpath, rawfile)); % load HCP rMEG data
    layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
else
    % If input file is raw eeg/meg file, this step is needed.
    % For HCP data, these steps are optional.
    %% Load MEG resting state data, Time windows of interset
    cfg                        = [];
    cfg.dataset                = fullfile(rmegpath, rawfile); 
    cfg.continuous             = 'yes';
    cfg.channel                = {upper(filetype)};
    
    data                       = ft_preprocessing(cfg);

    %% Redefine trial. Cut one time series data into short time series(trial)
    cfg                        = [];
    cfg.length                 = 2;
    cfg.overlap                = 0.5;
    
    data                       = ft_redefinetrial(cfg, data);
    
    %% Browse data.
    % You may need to check DC component and unnecessary channels to get 
    % rid of. And also, at the end of the recording, data may contain 0's.
    % So tou need to select trials.
    cfg                        = []; 
    [layoutpath layoutfile]    = fileparts(layoutfile); 
    cfg.layout                 = layoutfile;
    cfg.continuous             = 'no';
    cfg.viewmode               = 'vertical';
    ft_databrowser(cfg, data);
    
    %% Redefine data
    zerotrial                   = []; % Trials to get rid of.
    cfg                        = [];
    cfg.demean                 = 'yes';
    cfg.trials                 = setdiff(1:numel(data.trial), zerotrial);
    
    clear zerotrial;
    data                       = ft_preprocessing(cfg, data);
    %% Downsample
    cfg = [];
    cfg.resamplefs = 200; % Resample frequency set to 200 Hz
    cfg.detrend = 'yes';
    data = ft_resampledata(cfg, data);

    %% use ICA in order to identify cardiac and blink components
    cfg                 = [];
    cfg.method          = 'runica'; 
    cfg.runica.maxsteps = 50;
    cfg.randomseed      = 0; % this can be uncommented to match the data that has been stored on disk
    comp                = ft_componentanalysis(cfg, data);

    %% Reconstruct data from ica component
    cfg = [];
    cfg.component = [];% Specify badchannels;
    cfg.demean = 'no';
    cfg.updatesens = 'no';
    data = ft_rejectcomponent(cfg, comp);

end

%%     SPECTRUM ANALYSIS    %%
%% Calculate the Powersepctrum
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq =1;
cfg.keeptrials = 'no';
datapow = ft_freqanalysis(cfg, data);


%% plot powerspectrum
freqband = [2 4;5 7;8 12;15 30;30 60;60 90];
freqname = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma1', 'Gamme2'};

figure;
cfg = [];
cfg.layout = '4D248_helmet.mat';
for i =1:size(freqband, 1)
    cfg.xlim   = freqband(i,:);
    subplot(2,3,i);
    title([freqname{i} ' (' num2str(freqband(i,1)) '-' num2str(freqband(i,2)) 'Hz)']);
    ft_topoplotER(cfg, datapow);
    colorbar;
end

figure;
cfg = [];
cfg.channel = {'A129'}; % If multiple channel is set, then the function ft_singleplotER will calculate average power.
ft_singleplotER(cfg, datapow);

clear freqband freqname i;
%%    SOURCE RECONSTRUCTION OF SPECTRAL DATA    %%
%% Forward solution(Calculate leadfield)
headmodel = load(fullfile(megrootpath, anapath, headmodelfile));
headmodel = headmodel.headmodel;
sourcemodel = load(fullfile(megrootpath, anapath ,sourcemodelfile));
sourcemodel = sourcemodel.sourcemodel2d;

% % visualize the coregisteration of sensor, headmodel, and sourcemodel.
% figure;
% 
% % make the headmodel surface transparent
% ft_plot_vol(hdm, 'edgecolor', 'none'); alpha 0.4           
% ft_plot_mesh(ft_convert_units(sourcemodel, 'cm'),'vertexcolor',sourcemodel.sulc);
% ft_plot_sens(dataclean.grad);
% view([0 -90 0])

% Calculate Lead Field
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = headmodel;
cfg.channel = {'MEG'};
leadfield = ft_prepare_leadfield(cfg, data);

fid = fopen(fullfile(megrootpath,anapath,coordtransformfile), 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;

clear headmodelfile sourcemodelfile;
%% Frequency Analysis for Source Reconstruction
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foi = 10; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, data);

%% Source Reconstruction
cfg                   = [];
cfg.frequency         = spectdata.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, spectdata);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

%% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [1 8]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'yes';
ft_sourceplot(cfg, source);
view([-90 30]);
light;

%% compute sensor level single trial power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [8 12];                          
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, data);

%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));
%% compute the power spectrum for the median splitted data
cfg              = [];
cfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(cfg, datapow);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%% plot the topography of the difference along with the spectra
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [5 7];
figure; ft_topoplotER(cfg, powratio);

cfg         = [];
cfg.channel = {'A196'};
figure; ft_singleplotER(cfg, datapow_high, datapow_low);

%% compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;

cfg.trials = indlow; 
freq_low   = ft_freqanalysis(cfg, data);

cfg.trials = indhigh; 
freq_high  = ft_freqanalysis(cfg, data);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = spectdata.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, spectdata);

% use the precomputed filters 
cfg                   = [];
cfg.frequency         = spectdata.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.grid.filter       = source.avg.filter;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high, source_low);

% create a fancy mask
source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 

cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_ratio);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

%%
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foi = 10; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, data);

cfg                   = [];
cfg.frequency         = spectdata.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source                = ft_sourceanalysis(cfg, spectdata);

cfg                   = [];
cfg.method            = 'coh';
cfg.complex           = 'absimag';
source_connect        = ft_connectivityanalysis(cfg, source);

cfg = [];
cfg.method = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = 0.1;
network_full = ft_networkanalysis(cfg, source_connect);

cfg =[];
cfg.method = 'surface';
cfg.funparameter = 'degrees'
cfg.funcolormap = 'jet'
ft_sourceplot(cfg, network_full);
view([-150,30]);

sourcemodel_mni.tri = source.tri;
sourcemodel_mni.pos = source.pos*10;
sourcemodel_mni.pos(:,4) = 1;
sourcemodel_mni.pos = sourcemodel_mni.pos*transform.bti2spm';
sourcemodel_mni.pos = sourcemodel_mni.pos(:,1:3);
sourcemodel_mni.unit = 'mm';
sourcemodel_mni.coordsys = 'mni';

seedloc = struct('mni', zeros(4,3), 'label',repmat({''},4,1));
seedloc(1).mni = [-43,-76,35];
seedloc(1).label = 'LAG';
seedloc(1).net = 'DMN'
seedloc(2).mni = [51,-64,32];
seedloc(2).label = 'RAG';
seedloc(2).net = 'DMN'
seedloc(3).mni = [-3,-54,31];
seedloc(3).label = 'LPCC';
seedloc(3).net = 'DMN'
seedloc(4).mni = [-2,50.5,1.7];
seedloc(4).label = 'LMPFC';
seedloc(4).net = 'DMN'
seedloc(5).mni = [-25.3,-67.3,47.6];
seedloc(5).label = 'lpIPS';
seedloc(5).net = 'DAN'
seedloc(6).mni = [23.2,-69.4,48.6];
seedloc(6).label = 'rpIPS';
seedloc(6).net = 'DAN'
seedloc(7).mni = [-26.3,-11.8,52.7];
seedloc(7).label = 'lFEF';
seedloc(7).net = 'DAN'
seedloc(8).mni = [30.3,-12.8,52.6];
seedloc(8).label = 'rFEF';
seedloc(8).net = 'DAN'

for i=1:8
    distance_s2m(:,:,i) = zeros(8004,1);
    distance_s2m(:,:,i) = ( (sourcemodel_mni.pos(:,1) - seedloc(i).mni(1)).^2+...
                            (sourcemodel_mni.pos(:,2) - seedloc(i).mni(2)).^2+...
                            (sourcemodel_mni.pos(:,3) - seedloc(i).mni(3)).^2 ).^0.5;
    seedind(i).index  = find(distance_s2m(:,:,i) == min(distance_s2m(:,:,i)));
    seedind(i).dist   = distance_s2m(seedind(i).index, 1, i);
    seedloc(i).cohspctrm = source_connect.cohspctrm(seedind(i).index,:);
end

% source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 
source_plot = source;
% figure;
% imagesc(source_connect.cohspctrm)

%% start building the figure
h = figure('visible', 'on');
set(h, 'color', [1 1 1]);
set(h, 'renderer', 'opengl');
title('Seed based Functional Connectivity');

cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'coh';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [0 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';


surf     = [];
surf.pos = source_connect.pos;
surf.tri = source_connect.tri;




msk = zeros(size(seedloc(1).cohspctrm'));
msk(source_connect.inside) = 1;

% if hasroi && hasmsk
%   msk = roi .* msk;
%   opacmin = [];
%   opacmax = []; % has to be defined
% end

maskval = msk(:);
opacmin = 0;
opacmax = 1;
cfg.opacitymap = 'rampup';
cfg.opacitymap = alphamap(cfg.opacitymap);
alphamap(cfg.opacitymap);
%%
for i = 1:length(seedloc)
    val     = seedloc(i).cohspctrm';
    funmin = min(seedloc(i).cohspctrm(:));
    funmax = max(seedloc(i).cohspctrm(:));

    if isequal(cfg.funcolorlim, 'auto')
      if sign(funmin)>-1 && sign(funmax)>-1
        cfg.funcolorlim = 'zeromax';
      elseif sign(funmin)<1 && sign(funmax)<1
        cfg.funcolorlim = 'minzero';
      else
        cfg.funcolorlim = 'maxabs';
      end
    end
    if ischar(cfg.funcolorlim)
      % limits are given as string
      if isequal(cfg.funcolorlim, 'maxabs')
        fcolmin = -max(abs([funmin,funmax]));
        fcolmax =  max(abs([funmin,funmax]));
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'default'; end
      elseif isequal(cfg.funcolorlim, 'zeromax')
        fcolmin = 0;
        fcolmax = funmax;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'hot'; end
      elseif isequal(cfg.funcolorlim, 'minzero')
        fcolmin = funmin;
        fcolmax = 0;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'cool'; end
      else
        ft_error('do not understand cfg.funcolorlim');
      end
    else
      % limits are numeric
      fcolmin = cfg.funcolorlim(1);
      fcolmax = cfg.funcolorlim(2);
      % smart colormap
      if isequal(cfg.funcolormap, 'auto')
        if sign(fcolmin) == -1 && sign(fcolmax) == 1
          cfg.funcolormap = 'default';
        else
          if fcolmin < 0
            cfg.funcolormap = 'cool';
          else
            cfg.funcolormap = 'hot';
          end
        end
      end
    end
    clear funmin funmax

    figure;
    suptitle(['Seed based FC -' seedloc(i).net ', ' seedloc(i).label]);
    ft_plot_mesh(surf, 'edgecolor', 'none', 'facecolor', [], 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'facecolor', [], 'vertexcolor', val, 'facealpha', maskval, 'clim', [fcolmin fcolmax], 'alphalim', [opacmin opacmax], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');

    lighting gouraud
    view([-90 30]);
    light('style','infinite','position',[0 -200 200]);
    
end
%%
%% Source estimation with 'dics' and split trials into two groups
for k=1:1
    %% Frequent analysis
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq = 1;
    cfg.foilim = [8 12]; % set frequencies of interest as array
    spectdata = ft_freqanalysis(cfg, data);

    %% Source Reconstruction
    cfg                   = [];
    cfg.frequency         = [8,12];
    cfg.method            = 'dics';
    cfg.grid              = leadfield;
    cfg.headmodel         = headmodel;
    cfg.rawtrial          = 'yes';
    cfg.dics.lambda       = 0;
    cfg.dics.projectnoise  = 'yes';
    source = ft_sourceanalysis(cfg, spectdata);

    %% Find ROI
    use_my_library('ft',0);
    use_my_library('spm',1);
    btiROI = spm_eeg_inv_transform_points(transform.spm2bti,[0,-58,0]); %% pcc
    use_my_library('spm',0);
    use_my_library('ft',1);
    tempDistance = sourcemodel.pos - repmat(btiROI,size(sourcemodel.pos,1),1);
    tempDistance = sum(tempDistance.^2,2);
    roiInd = find(min(tempDistance)==tempDistance);
    clear tempDistance btiROI
    %% Split trials into two group
    tmp = [];
    for i=1:length(source.trial)
        tmp(end+1) = source.trial(i).pow(roiInd);
    end
    indlow  = find(tmp(:)<=median(tmp(:)));
    indhigh = find(tmp(:)>=median(tmp(:)));
    %%
    cfg            = [];
    cfg.method     = 'mtmfft';
    cfg.output     = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq  = 1;
    cfg.foi        = 10;
    cfg.trials = indlow; 
    freq_low   = ft_freqanalysis(cfg, data);
    cfg.trials = indhigh; 
    freq_high  = ft_freqanalysis(cfg, data);
    
    cfg                   = [];
    cfg.frequency         = [8 12];
    cfg.method            = 'dics';
    cfg.grid              = leadfield;
    cfg.headmodel         = headmodel;
    cfg.keeptrials        = 'yes';
    cfg.dics.lambda        = '10%';
    cfg.dics.projectnoise  = 'yes';
    cfg.dics.keepfilter    = 'yes';
    source = ft_sourceanalysis(cfg, spectdata);

    % use the precomputed filters 
    cfg                   = [];
    cfg.frequency         = [8,12];
    cfg.method            = 'dics';
    cfg.grid              = leadfield;
    cfg.grid.filter       = source.avg.filter;
    cfg.headmodel         = headmodel;
    cfg.keeptrials        = 'yes';
    cfg.dics.lambda        = '10%';
    cfg.dics.projectnoise  = 'yes';
    source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
    source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

    cfg           = [];
    cfg.operation = 'log10(x1)-log10(x2)';
    cfg.parameter = 'pow';
    source_ratio  = ft_math(cfg, source_high, source_low);

    % create a fancy mask
    source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 

    cfg = [];
    cfg.method        = 'surface';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = 'mask';
    cfg.funcolorlim   = [-.3 .3];
    cfg.funcolormap   = 'jet';
    cfg.colorbar      = 'no';
    ft_sourceplot(cfg, source_ratio);
    view([-90 30]);
    light('style','infinite','position',[0 -200 200]);
    %%
    cfg            = [];
    cfg.method     = 'mtmfft';
    cfg.output     = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq  = 1;
    cfg.foi        = 10;
    cfg.trials = indlow; 
    freq_low   = ft_freqanalysis(cfg, data);
    cfg.trials = indhigh; 
    freq_high  = ft_freqanalysis(cfg, data);
    
    cfg                   = [];
    cfg.frequency         = [8 12];
    cfg.method            = 'pcc';
    cfg.grid              = leadfield;
    cfg.headmodel         = headmodel;
    cfg.keeptrials        = 'yes';
    cfg.pcc.lambda        = 0;
    cfg.pcc.projectnoise  = 'yes';
    cfg.pcc.keepfilter    = 'yes';
    cfg.pcc.fixedori      = 'yes';
    source = ft_sourceanalysis(cfg, spectdata);

    % use the precomputed filters 
    cfg                   = [];
    cfg.frequency         = [8,12];
    cfg.method            = 'pcc';
    cfg.grid              = leadfield;
    cfg.grid.filter       = source.avg.filter;
    cfg.headmodel         = headmodel;
    cfg.keeptrials        = 'yes';
    cfg.pcc.lambda        = 0;
    cfg.pcc.projectnoise  = 'yes';
    source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
    source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

    cfg           = [];
    cfg.operation = 'log10(x1)-log10(x2)';
    cfg.parameter = 'pow';
    source_ratio  = ft_math(cfg, source_high, source_low);

    % create a fancy mask
    source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 

    cfg = [];
    cfg.method        = 'surface';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = 'mask';
    cfg.funcolorlim   = [-.3 .3];
    cfg.funcolormap   = 'jet';
    cfg.colorbar      = 'no';
    ft_sourceplot(cfg, source_ratio);
    view([-90 30]);
    light('style','infinite','position',[0 -200 200]);
end