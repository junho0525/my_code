function mnet_ft_sourceplot(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mnet source power visualize using ft_freqanalysis, ft_sourceanalysis,
%   and ft_sourceplot
% 
% Perform multitaper fourier transform useing ft_freqanalysis and then
% visualize sensor space power spectrum
%
% USAGE
%   mnet_ft_topoplot(data, band, anafile, datalabels, plotsize)
%     data : fieldtrip struct or cell array of fieldtrip struct.
%     band : 'alpha', 'beta', 'delta','theta','gamma1','gamma2'
%          and any frequency band or cell array of frequency bands of
%          interest.
%     datalabels(Optional) : if you want to specify data name on each
%                         plot, make datalabels as the same size of data.
%     plotsize (Optional) : size of row and column that will be used for
%                         subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.09.18 11:47 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
use_my_library init
use_my_library('ft',1);
if numel(varargin)>5
    
end
data = varargin{1};
band = varargin{2};
anafile = varargin{3};
if iscell(data), ismultisubj=1; else, ismultisubj = 0; subjnum = 1; end
if ismultisubj, subjnum = numel(data); end
if ismultisubj && length(anafile)==1
    warning('Use single anatomy data for all MEG data');
    error('This is not implemented yet');
elseif ismultisubj && length(anafile)>1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
if numel(varargin)<4
    if ismultisubj
        datalabels = strcat('Data\_',strsplit(num2str(1:subjnum),' '));
    else
        datalabels = 'Data';
    end
else
    datalabels = varargin{4};
    if numel(data)~=numel(datalabels), error('data and label do not match in numbers.'); end
end
if iscell(band), ismultiband=1; else, ismultiband = 0; end
if ismultiband 
    bandnum = numel(band);
    for i = 1:bandnum
        if ischar(band{i})
            band{i} = lower(band{i});
            switch band{i}
                case 'delta'
                    bandfreq{i} = [2 4];
                case 'theta'
                    bandfreq{i} = [5 7];
                case 'alpha'
                    bandfreq{i} = [8 12];
                case 'beta'
                    bandfreq{i} = [15 30];
                case 'gamma1'
                    bandfreq{i} = [30 60];
                case 'gamma2'
                    bandfreq{i} = [60 90];
                case 'all'
                    bandfreq{i} = [0 data{i}.fsample/2];
            end
        end
    end
else
    bandnum = 1;
    if ischar(band)
        band = lower(band);
        switch band
            case 'delta'
                bandfreq= [2 4];
            case 'theta'
                bandfreq = [5 7];
            case 'alpha'
                bandfreq = [8 12];
            case 'beta'
                bandfreq = [15 30];
            case 'gamma1'
                bandfreq = [30 60];
            case 'gamma2'
                bandfreq = [60 90];
            case 'all'
                bandfreq{i} = [0 data{i}.fsample/2];
        end
    end
end
if numel(varargin)>3
    subplotsize = varargin{4};
else
    subplotsize = [];
end
%% Power anlaysis of sensor space
datapow = [];
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq =1;
cfg.keeptrials = 'no';
for i = 1:subjnum
    if isfield(data{i}, 'data')
        datapow{i} = ft_freqanalysis(cfg, data{i}.data);
    else
        datapow{i} = ft_freqanalysis(cfg, data{i});
    end
end

%% Plot sensor space power spectrum
monitorpos = get(0,'MonitorPosition');
if isempty(subplotsize)
    switch bandnum
        case 1
            plotnum = 12;
            subplotsize = [3,4];
        case 2
            plotnum = 12;
            subplotsize = [3,4];
        otherwise
            plotnum = 12;
            subplotsize = [3, bandnum];
    end
    for i = 1:numel(datapow)
        if mod((i-1),plotnum)==0
            figure('outerposition',monitorpos(1,:));
        end
        plotindex = mod((i-1),plotnum)+1;
        subplot(subplotsize(1), subplotsize(2), plotindex);
        
        cfg = [];
        cfg.layout = '4D248_helmet.mat';
        cfg.xlim   = bandfreq;
        title([datalabels{i} '\_' band]);
        ft_topoplotER(cfg, datapow{i});
        colorbar;
    end
else
    for i = 1:numel(datapow)
        plotnum = subplotsize(1)*subplotsize(2);
        if mod((i-1),plotnum)==0
            figure('outerposition',monitorpos(1,:));
        end
        plotindex = mod((i-1),plotnum)+1;
        subplot(subplotsize(1), subplotsize(2), plotindex);
        
        cfg = [];
        cfg.layout = '4D248_helmet.mat';
        cfg.xlim   = bandfreq;
        title([datalabels{i} '\_' band]);
        ft_topoplotER(cfg, datapow{i});
        colorbar;
    end
end