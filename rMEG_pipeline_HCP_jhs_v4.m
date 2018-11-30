function D = rMEG_pipeline_HCP_jhs_v4(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP Task MEG Source Reconstruction Pipeline Using FieldTrip             %
%                                                                         %
% This pipeline code is for rest HCP dataset                              %
%                                                                         %
% D = rMEG_pipeline_HCP_jhs_v4(config)                                    %
%    Output      D     - spm meeg object                                  %
%    Input                                                                %
%       config.subjectID     - HCP dataset subject id                     %
%       config.tasktype      - 'Restin' (Default)                         %
%       config.subjectPath   - HCP Single Subject Data Path               %
%       config.savePath      - Save path for output and intermediate files%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.04.02 20:22 - By Junho Son                                     %
%         - Make this code based on the code rMEG_pipeline_HCP_jhs_v2.0.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end

% These fields is optional. If these fields are not specified, use default.
if ~isfield(config, 'taskType'),config.taskType = 'Restin';end % default : 'Restin'
if ~isfield(config, 'sourceLocalization'), config.sourceLocalization = 1; end % default : 1 
if ~isfield(config, 'showForwardResult'), config.showForwardResult = 0; end % default : 0
if ~isfield(config, 'showMEGSignal'), config.showMEGSignal = 0; end % default : 0
if ~isfield(config, 'showAlphaPowRatio'), config.showAlphaPowRatio = 0; end % default : 0
if ~isfield(config, 'showConnectionDensityMap'), config.showConnectionDensityMap = 0; end % default : 0
if ~isfield(config, 'showFullFC'), config.showFullFC = 0; end % default : 0
if ~isfield(config, 'showSeedBasedFC'), config.showSeedBasedFC = 1; end % default : 1
if ~isfield(config, 'saveFig'), config.saveFig = 0; end % default : 0
if ~isfield(config, 'showSensPower'), config.showSensPower = 0; end % default : 0
%% Set Path
switch config.taskType
    case 'Restin'
        rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
        rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
        icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-Restin_icaclass_vs.mat']);
    otherwise
        Error();
end
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
coordtransformfile = [config.subjectID '_MEG_anatomy_transform.txt'];
% Set Save Path
if config.saveFig
    if ~isdir(fullfile(config.savePath, 'figure')), mkdir(fullfile(config.savePath, 'figure'));end
    if ~isdir(fullfile(config.savePath, 'figure',[config.subjectID '_' num2str(config.sessionID)])), mkdir(fullfile(config.savePath, 'figure',[config.subjectID '_' num2str(config.sessionID)]));end
    savePath = fullfile(config.savePath, 'figure',[config.subjectID '_' num2str(config.sessionID)]);
end
%%     PREPROCESSING     %%
%% Load MEG data
load(fullfile(rMEGPath,rMEGRawData)); % load HCP rMEG data
layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
try
    load(fullfile(config.savePath, ['LF_' config.subjectID '_' num2str(config.sessionID) '_Restin']));
    isLeadField = 1;
catch
    isLeadField = 0;
end
try
    load(fullfile(config.savePath, ['SourceConn_' config.subjectID '_' num2str(config.sessionID) '_Restin']));
    isSourceConn = 1;
catch
    isSourceConn = 0;
end
%% Browse data.
if config.showMEGSignal
    cfg                    = []; 
    cfg.layout             = layoutfile;
    cfg.continuous         = 'no';
    cfg.viewmode           = 'vertical';
    ft_databrowser(cfg, data);
end
%%     SPECTRUM ANALYSIS    %%
if config.showSensPower
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
        if config.saveFig
            if ~isdir(fullfile(savePath, 'SensPower')),mkdir(fullfile(savePath, 'SensPower'));end
            saveas(gcf ,fullfile(savePath, 'SensPower',['FT_' config.subjectID '_' num2str(config.sessionID) '_sensPow.tif']));
        end
    end
    clear freqband freqname i datapow;
end
%%    CALCULATE LEADFIELD    %%
%% Forward solution(Calculate leadfield)
headmodel = load(headmodelfile);
headmodel = headmodel.headmodel;
sourcemodel = load(sourcemodelfile);
sourcemodel = sourcemodel.sourcemodel2d;
if ~isLeadField
    cfg = [];
    cfg.grid = sourcemodel;
    cfg.headmodel = headmodel;
    cfg.channel = {'MEG'};
    leadfield = ft_prepare_leadfield(cfg, data);
    save(fullfile(config.savePath, ['LF_' config.subjectID '_' num2str(config.sessionID) '_Restin']),'leadfield');
end
fid = fopen(fullfile(anaPath,coordtransformfile), 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);

clear fid strFid;
clear headmodelfile sourcemodelfile;


%%    SPECTRAL ANALYSIS - ALPHA POWER RATIO    %%
if config.showAlphaPowRatio
     %% Frequency Analysis for Source Reconstruction
%     cfg = [];
%     cfg.method = 'mtmfft';
%     cfg.output = 'fourier';
%     cfg.keeptrials = 'yes';
%     cfg.tapsmofrq = 1;
%     cfg.foi = 10; % set frequencies of interest as array
%     spectdata = ft_freqanalysis(cfg, data);
%     %% Source Reconstruction
%     cfg                   = [];
%     cfg.frequency         = spectdata.freq;
%     cfg.method            = 'pcc';
%     cfg.grid              = leadfield;
%     cfg.headmodel         = headmodel;
%     cfg.keeptrials        = 'yes';
%     cfg.pcc.lambda        = '10%';
%     cfg.pcc.projectnoise  = 'yes';
%     cfg.pcc.fixedori      = 'yes';
%     source = ft_sourceanalysis(cfg, spectdata);
%     source = ft_sourcedescriptives([], source); % to get the neural-activity-index
%      %% plot the neural activity index (power/noise)
%     cfg               = [];
%     cfg.method        = 'surface';
%     cfg.funparameter  = 'nai';
%     cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = 'auto'; % [0.0 8]
%     cfg.opacitylim    = [3 8]; 
%     cfg.opacitymap    = 'rampup';  
%     cfg.funcolormap   = 'jet';
%     cfg.colorbar      = 'yes';
%     ft_sourceplot(cfg, source);
%     view([-90 30]);
%     light;
%     clear source
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
    clear tmp freqind chanind
    %% compute the power spectrum for the median splitted data
    cfg              = [];
    cfg.trials       = indlow; 
    datapow_low      = ft_freqdescriptives(cfg, datapow);

    cfg.trials       = indhigh; 
    datapow_high     = ft_freqdescriptives(cfg, datapow);
    %% compute the difference between high and low
%     cfg = [];
%     cfg.parameter = 'powspctrm';
%     cfg.operation = 'divide';
%     powratio      = ft_math(cfg, datapow_high, datapow_low);

    %% plot the topography of the difference along with the spectra
%     cfg        = [];
%     cfg.layout = '4D248_helmet.mat';
%     cfg.xlim   = [9.9 10.1];
%     figure; ft_topoplotER(cfg, powratio);
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
    clear indlow indhigh
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

    %% create a fancy mask
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
    if config.saveFig
        if ~isdir(fullfile(savePath, 'AlphaPowRatio')),mkdir(fullfile(savePath, 'AlphaPowRatio'));end
        saveas(gcf, fullfile(savePath, 'AlphaPowRatio', ['FT_' config.subjectID '_' num2str(config.sessionID) '_AlphaRatio.tif'] ))
    end
end
%%    SPECTRAL ANALYSIS - SPECTRAL FC ANALYSIS    %%
if ~isSourceConn
    %% Frequency Analysis 
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
    cfg.pcc.keepfilter    = 'yes';
    cfg.pcc.fixedori      = 'yes';
    source                = ft_sourceanalysis(cfg, spectdata);
    %% Calculate Functional Connectivity using Coherence
    cfg                   = [];
    cfg.method            = 'coh';
    cfg.complex           = 'absimag';
    source_connect        = ft_connectivityanalysis(cfg, source);
    save(fullfile(config.savePath, ['SourceConn_' config.subjectID '_' num2str(config.sessionID) '_Restin']),'source_connect');
end
%% Network Analysis - show rich club nodes
if config.showConnectionDensityMap
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
    if config.saveFig
        if ~isdir(fullfile(savePath, 'ConnDensityMap')),mkdir(fullfile(savePath, 'ConnDensityMap'));end
        saveas(gcf, fullfile(savePath, 'ConnDensityMap', ['FT_' config.subjectID '_' num2str(config.sessionID) '_ConnDensity.png'] ))
    end
end
%% Show Full Functional Connectivity Matrix
if config.showFullFC
    figure;
    imagesc(source_connect.cohspctrm)
    if config.saveFig
        if ~isdir(fullfile(savePath, 'FullFC')),mkdir(fullfile(savePath, 'FullFC'));end
        saveas(gcf, fullfile(savePath, 'FullFC', ['FT_' config.subjectID '_' num2str(config.sessionID) '_FullFC_coh.png'] ))
    end
end
%% 
if config.showSeedBasedFC
    sourcemodel_mni.tri = source_connect.tri;
    sourcemodel_mni.pos = source_connect.pos*10;
    sourcemodel_mni.pos(:,4) = 1;
    sourcemodel_mni.pos = sourcemodel_mni.pos*transform.bti2spm';
    sourcemodel_mni.pos = sourcemodel_mni.pos(:,1:3);
    sourcemodel_mni.unit = 'mm';
    sourcemodel_mni.coordsys = 'mni';

    seedloc = struct('mni', zeros(4,3), 'label',repmat({''},4,1));
    seedloc(1).mni = [-43,-76,35];
    seedloc(1).label = 'LAG';
    seedloc(1).net = 'DMN';
    seedloc(2).mni = [51,-64,32];
    seedloc(2).label = 'RAG';
    seedloc(2).net = 'DMN';
    seedloc(3).mni = [-3,-54,31];
    seedloc(3).label = 'LPCC';
    seedloc(3).net = 'DMN';
    seedloc(4).mni = [-2,50.5,1.7];
    seedloc(4).label = 'LMPFC';
    seedloc(4).net = 'DMN';
    seedloc(5).mni = [-25.3,-67.3,47.6];
    seedloc(5).label = 'lpIPS';
    seedloc(5).net = 'DAN';
    seedloc(6).mni = [23.2,-69.4,48.6];
    seedloc(6).label = 'rpIPS';
    seedloc(6).net = 'DAN'
    seedloc(7).mni = [-26.3,-11.8,52.7];
    seedloc(7).label = 'lFEF';
    seedloc(7).net = 'DAN';
    seedloc(8).mni = [30.3,-12.8,52.6];
    seedloc(8).label = 'rFEF';
    seedloc(8).net = 'DAN';
    
    % Calculate Euclidean distance from seed to all other points and find
    % the nearest point for new seed on the surface
    for i=1:8
        distance_s2m(:,:,i) = zeros(8004,1);
        distance_s2m(:,:,i) = ( (sourcemodel_mni.pos(:,1) - seedloc(i).mni(1)).^2+...
                                (sourcemodel_mni.pos(:,2) - seedloc(i).mni(2)).^2+...
                                (sourcemodel_mni.pos(:,3) - seedloc(i).mni(3)).^2 ).^0.5;
        seedind(i).index  = find(distance_s2m(:,:,i) == min(distance_s2m(:,:,i)));
        seedind(i).dist   = distance_s2m(seedind(i).index, 1, i);
        seedloc(i).cohspctrm = source_connect.cohspctrm(seedind(i).index,:);
    end
    clear distance_s2m sourcemodel_mni seedind
    %% Render FC on Cortical Surface
    %% This code is Originally from ft_sourceplot.m
    %  configure for plot
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
    %%%% This is for ROI based Analysis
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
    %% Plot FC for Each Seeds
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
        h = figure;
        set(h, 'color', [1 1 1]);
        set(h, 'renderer', 'opengl');
        suptitle(['Seed based FC -' seedloc(i).net ', ' seedloc(i).label]);
        ft_plot_mesh(surf, 'edgecolor', 'none', 'facecolor', [], 'vertexcolor', 'curv');
        ft_plot_mesh(surf, 'edgecolor', 'none', 'facecolor', [], 'vertexcolor', val, 'facealpha', maskval, 'clim', [fcolmin fcolmax], 'alphalim', [opacmin opacmax], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
        lighting gouraud
        light('style','infinite','position',[0 -200 200]);
        light('style','infinite','position',[0 200 -200]);
        drawnow
        if config.saveFig
            if ~isdir(fullfile(savePath, 'SeedFC')), mkdir(fullfile(savePath, 'SeedFC'));end
            if ~isdir(fullfile(savePath, 'SeedFC', 'DAN')), mkdir(fullfile(savePath, 'SeedFC', 'DAN'));end
            if ~isdir(fullfile(savePath, 'SeedFC', 'DMN')), mkdir(fullfile(savePath, 'SeedFC', 'DMN'));end
            viewList = [];
            viewList = {[1,0,0],[-1,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,-1,0],[0,1,0]};
            figureLabel = {'A', 'P','S','I','L','LM','R','RM'};
            for j =1:8
                if j > 4 && j <= 6
                    set(gca, 'YLim', [0,10], 'XLim', [-10 10], 'ZLim',[-10 10]);
                elseif j>6 && j<= 8
                    set(gca, 'YLim', [-10,0],'XLim', [-10 10], 'ZLim',[-10 10]);
                end
                view(viewList{j});
                if j ==3 || j==4
                    camroll(90); % upside of fig is anterior
                end
                drawnow
                t = timer;
                t.StartDelay = 0.1;
                t.TimerFcn = @(~,~) print(fullfile(savePath, 'SeedFC', seedloc(i).net,  ['FT_HCP_MEG_' config.subjectID '_' num2str(config.sessionID) '_SeedFC_' seedloc(i).label '_' figureLabel{j}] ), '-dtiff', '-noui');
                start(t)
                wait(t)
                delete(t)
                clear t
            end
            close(h);
            
        end
    end
    clear surf msk seedloc maskval opacmin opacmax val i fcolmin fcolmax cfg
end
D = [];