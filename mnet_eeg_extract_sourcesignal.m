function sourceSignal = mnet_eeg_extract_sourcesignal(meegFile       , ...
                                             processMode,options,saveInfo)
%   MNET_EEG_EXTRACT_SORUCESIGNAL extract source space timeseries using
%   "Headmodel" and "Sourcemodel". If you have own headmodel or sourcemodel
%   please read process mode.
%
%   ======================================================================
%   # Input
%   ======================================================================
%      -------------------------------------------------------------------
%      EEG file            [E]   : EEG filename.
%      -------------------------------------------------------------------
%      Process mode        [O]   : Process on/off setting.
%      -------------------------------------------------------------------
%       * headmodel              | Make headmodel or Use individual model.
%       * sourcemodel            | Make source model or Use individual
%                                  model.
%       * downsample             | Downsampling data.
%       * filtering              | Bandpass filtering.
%      -------------------------------------------------------------------
%      Options             [O]   : Optional values.
%                                  If empty, default value will assign.
%      -------------------------------------------------------------------
%       * headmodel              | Headmodel file or Individual sMRI.
%       * sourcemode             | Sourcemodel file                 .
%       * samplingrate           | Down sampling rate.
%       * bandpass               | Filtering band.
%      -------------------------------------------------------------------
%   ======================================================================
%   # Output
%   ======================================================================
%      -------------------------------------------------------------------
%      Source signal          : Source-space signal
%      -------------------------------------------------------------------
%   ======================================================================
%   Examples
%   ======================================================================
%   ______________________________________________________________________
%   $ 2018.10.14 : Initial Release.

%%  0. Varargin Check
%   ======================================================================
%   0.1 If no variabel
%   ----------------------------------------------------------------------
if  nargin < 1, meegFile = spm_select; end
if  nargin < 2, processMode = struct; end
if  nargin < 3, options     = struct; end
if  nargin < 4, saveInfo    = struct; end
%   ----------------------------------------------------------------------
%   0.2 Process Mode
%       : Process On/Off
%   ----------------------------------------------------------------------
if  ~(isfield(processMode,'headmodel'))  , processMode.headmodel   = 0;end
if  ~(isfield(processMode,'sourcemodel')), processMode.sourcemodel = 0;end
if  ~(isfield(processMode,'downsample')) , processMode.downsample  = 1;end
if  ~(isfield(processMode,'filtering'))  , processMode.filtering   = 1;end
%   ----------------------------------------------------------------------
if  isfield(options,'headmodel')
    processMode.headmodel = 1; 
else
    processMode.headmodel = 0;
end
if  isfield(options,'sourcemodel')
    processMode.sourcemodel = 0;
end
%   ----------------------------------------------------------------------
%   0.3 Options
%   ----------------------------------------------------------------------
if  ~(isfield(options,'bandpass'))       , options.bandpass =[0.1 100];end
if  ~(isfield(options,'samplingrate'))   , options.samplingrate = 250 ;end
%   ----------------------------------------------------------------------
%   0.4 Save Info
%   ----------------------------------------------------------------------
if  ~(isfield(saveInfo,'resultfile'))    , saveInfo.resultfile     ='';end
%   ----------------------------------------------------------------------
%   0.5 Path setting
%   ----------------------------------------------------------------------
    check = which('pop_loadset');
if  isempty(check)
    error('Addpath EEGLab subfolder');
end
    check = which('spm_eeg_load');
if  isempty(check)
    error('Addpath SPM EEG folder');
end
    check = which('ft_prepare_leadfield');
if  isempty(check)
    error('Addpath FieldTrip folder');
end
%   ----------------------------------------------------------------------
%   0.6 Check type 
%   ----------------------------------------------------------------------
if  isempty(meegFile)
    error('\t - EEG file is empty');
end
    [~,~,ext] = fileparts(meegFile);
%   ----------------------------------------------------------------------
%   - EEGLab format    
%   ----------------------------------------------------------------------
%   * Load EEGLab file ( pop_loaset, set file )
if  strcmpi(ext,'.fdt')
    meegFile = spm_file(meegFile,'ext','set');
end
if  strcmpi(ext,'.dat')
    meegFile = spm_file(meegFile,'ext','mat');
end
if  strcmpi(ext,'.set') 
    EEG  = pop_loadset(meegFile);
    data = eeglab2fieldtrip(EEG,'preprocessing');

%   ----------------------------------------------------------------------
%   - SPM format
%   ----------------------------------------------------------------------
elseif  strcmpi(ext,'.mat')
    if  exist(spm_file(meegFile,'ext','dat'),'file')
        D    = spm_eeg_load(meegFile);
    %   * Convert to FieldTrip format
        data         = D.ftraw;
        data.fsample = D.fsample;
    else
%   ----------------------------------------------------------------------
%   - FieldTrip format
%   ----------------------------------------------------------------------
        load(meegFile);
    end
else
    error('Unsupproted format');
end
%   ----------------------------------------------------------------------
%   - M/EEG Check
%   ----------------------------------------------------------------------
if  isfield(data,'elec')
    options.type = 'eeg';
end
if  isfield(data,'grad')
    options.type = 'meg';
end

%%  1. Preprocess before Source reconstruction
%   ======================================================================
%   1.1 Downsampling
%   ----------------------------------------------------------------------
if  processMode.downsample && (data.fsample > options.samplingrate)
    fprintf(1,'\t - Downsampling %d to %d\n',data.fsample,     ...
                                             options.samplingrate);
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'yes';
    cfg.demean     = 'yes';
    data           = ft_resampledata(cfg,data);
end
%   ----------------------------------------------------------------------
%   1.2 Filtering
%       @ Must satisfy this condition : (Freqmax * 2 < fsample) 
%   ----------------------------------------------------------------------
if  processMode.filtering
    fprintf(1,'\t - Bandpass filtering 0.1 to 100 Hz\n');
    cfg            = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [0.1 100];
    data           = ft_preprocessing(cfg,data);
end
%%  2. Load or Make Headmodel and Sourcemodel
%   ======================================================================
%   2.1 Head model
%   ----------------------------------------------------------------------
if  processMode.headmodel
%   + To be update
    
else
    ftPath = fileparts(check);
if  strcmpi(options.type','meg')
    options.headmodel = fullfile(ftPath,'template','headmodel', ...
                                    'standard_singleshell.mat'      );
else
    options.headmodel = fullfile(ftPath,'template','headmodel', ...
                                            'standard_bem.mat'      );
end
end
if  ~(exist(options.headmodel,'file'))
    error('Headmodel file is not exist.');
else
    load(options.headmodel);
    headmodel = vol;
end
%   ----------------------------------------------------------------------
%   2.2 Source model
%   ----------------------------------------------------------------------
if  processMode.sourcemodel
    
else
    ftPath = fileparts(check);
    options.sourcemodel = fullfile(ftPath,'template','sourcemodel', ...
                                                 'cortex_8196.surf.gii');
end
if  ~(exist(options.sourcemodel))
    error('Sourcemodel file is not exist');
else
    sourcemodel = ft_read_headshape(options.sourcemodel);
end
%%  3. Construct Leadfield and Source reconstruction
%   ======================================================================
%   3.1 Construct leadfield matrix
%   ----------------------------------------------------------------------
    fprintf(1,'\t - Construct Leadfield matrix\n');
    cfg            = [];
    cfg.headmodel  = headmodel;
    cfg.grid       = sourcemodel;
if  strcmpi(options.type,'meg')
    cfg.grad       = data.grad;
else
    cfg.elec       = data.elec;
end
    cfg.channel    = data.label;
    leadfield      = ft_prepare_leadfield(cfg);
%   ----------------------------------------------------------------------
%   3.2 Source Reconstruction (using Beamformer LCMV)
%   ----------------------------------------------------------------------
    fprintf(1,'\t - Source reconstruction\n');
    cfg            = [];
    cfg.grid       = leadfield;
    cfg.headmodel  = headmodel;
    cfg.method     = 'lcmv';
    
    cfg.lcmv.labmda       = '10%';
    cfg.lcmv.fixedori     = 'yes';
    cfg.lcmv.projectmoise = 'yes';
    
    cfg.rawtrial   = 'yes';      %@ Compute source signals for every trials
    source         = ft_sourceanalysis(cfg, data);

for i = 1:length(source.trial)
    for j = 1:length(source.trial(i).mom)
        if isempty(source.trial(i).mom{j})
            sourceSignal(j,:,i) = zeros(1,length(source.time));
        else
            sourceSignal(j,:,i)= source.trial(i).mom{j};
        end
    end
end
    fprintf(1,'\t * Source reconstruction done.\n');
%%  4. Save
%   ======================================================================
%   + To be update
%   ___________________________End of Function____________________________

    