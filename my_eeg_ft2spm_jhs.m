function D = my_eeg_ft2spm_jhs(filename,savePath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert fieldtrip task data to SPM MEEG object                          %
%                                                                         %
% This code for task HCP dataset                                          %
%                                                                         %
% D = my_eeg_ft2spm_sjh(ftdata, type, filename, savePath)                 %
%    Output      D     - spm meeg object                                  %
%    Input                                                                %
%        
%                                                                         %
% This code is modified from code spm_eeg_ft2spm.m (SPM12)                %
% Vladimir Litvak                                                         %
% $Id: spm_eeg_ft2spm.m 5438 2013-04-24 10:38:47Z vladimir $              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
% -- 2018.02.20 16:31 -- By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
%%
%  Input  filename
%  Result taskType, taskSubType, sessionNum, trialinfo, pathname, fname

isTF = 0;

% If raw format ; check data\
if isempty(filename)
    error('Invalid filename');
end
[pathname, fname] = fileparts(filename);
filenameparts = strsplit(fname, '_');
taskType = filenameparts{3};
taskType = strsplit(taskType,'-');
taskType = taskType{2};
clear filenameparts;

%% Load HCP Data and Extract Some Informations
ftdata = load(filename);

if isfield(ftdata, 'data')
    ftdata = ftdata.data;
end

if iscell(ftdata.time)
    Ntrials=length(ftdata.trial);
    Nchannels  = size(ftdata.trial{1},1);
    Nsamples  = size(ftdata.trial{1},2);

    data = zeros(Nchannels, Nsamples, Ntrials);
    % Because SPM can only handle data with same trial length and same 
    % trial borders, some task should redefine some trials.
    % Eg) Motort task fixation time window is [0s, 1.2s] but others are
    % [-1.2s, 1.2s]. So, We use zero-padding for non-exist fixation data.
    switch taskType
        case 'Wrkmem'
            for n=1:Ntrials
                data(:,:,n) = ftdata.trial{n};
            end
        case 'Motort'
            for n=1:Ntrials
                if ftdata.trialinfo(n,2) == 6
                    data(:,(Nsamples-length(ftdata.trial{n})+1):Nsamples,n) = ftdata.trial{n};
                else
                    data(:,:,n) = ftdata.trial{n};
                end
            end
        case 'StoryM'
        otherwise
    end
    
    % Convert data unit from T to fT, multiply 10e15
    data = data*1.0e15;
    
    ftdata.time  = ftdata.time{1};
else
    Nchannels  = numel(ftdata.label);
    Nsamples   = length(ftdata.time);
    
    rptind=strmatch('rpt', tokenize(ftdata.dimord, '_'));
    if isempty(rptind)
        rptind=strmatch('subj', tokenize(ftdata.dimord, '_'));
    end

    timeind=strmatch('time', tokenize(ftdata.dimord, '_'));
    chanind=strmatch('chan', tokenize(ftdata.dimord, '_'));

    if any(ismember({'trial', 'individual', 'avg'}, fieldnames(ftdata) )) % timelockanalysis
        if ~isempty(rptind)
            if isfield(ftdata, 'trial')
                Ntrials = size(ftdata.trial, rptind);
                data =permute(ftdata.trial, [chanind, timeind, rptind]);
            else
                Ntrials = size(ftdata.individual, rptind);
                data =permute(ftdata.individual, [chanind, timeind, rptind]);
            end
        else
            Ntrials = 1;
            data =permute(ftdata.avg, [chanind, timeind]);
        end        
    elseif isfield(ftdata, 'powspctrm')
        isTF = 1;
        Nfrequencies = numel(ftdata.freq);
        freqind = strmatch('freq', tokenize(ftdata.dimord, '_'));
        if ~isempty(rptind)
            Ntrials = size(ftdata.powspctrm, rptind);
            data = permute(ftdata.powspctrm, [chanind, freqind, timeind, rptind]);
        else
            Ntrials = 1;
            data = permute(ftdata.powspctrm, [chanind, freqind, timeind]);
        end
    end
end

%%
%--------- Start making the header

D = [];

% sampling rate in Hz
if isfield(ftdata, 'fsample')
    D.Fsample = ftdata.fsample;
else
    D.Fsample = 1./mean(diff(ftdata.time));
end

D.timeOnset = ftdata.time(1);

% Number of time bins in peri-stimulus time
D.Nsamples = Nsamples;

% Specify channels
% Set channel labels and units
% unit was converted from T to fT in Line 40. Multiplied by 10e15
D.channels = struct('label', ftdata.label,'units',repmat({'fT'},Nchannels,1));

D.trials = repmat(struct('label', {'Undefined'}), 1, Ntrials);
switch taskType
    case 'Wrkmem'
        % trialinfo(:,4) - 1: face 2: tools 0: fixation
        % trialinfo(:,5) - 1: 0-back 2: 2-back NaN: fixation
        % We define conditions only 0 back and 2 back for now.

        for n = 1:Ntrials
            switch ftdata.trialinfo(n,5)
                case 1
                    D.trials(n).label = '0-back';
                case 2
                    D.trials(n).label = '2-back';
                otherwise
                    D.trials(n).label = 'Fixation';
            end
        end
    case 'Motort'
        % ftdata.trialinfo(:,2) - 1: Left Hand 2: Left Foot 4: Right Hand
        % 5: Right Foot 6: Fixation
        for n = 1:Ntrials
            switch ftdata.trialinfo(n,2)
                case 1
                    D.trials(n).label = 'Left Hand';
                case 2
                    D.trials(n).label = 'Left Foot';
                case 4
                    D.trials(n).label = 'Right Hand';
                case 5
                    D.trials(n).label = 'Right Foot';
                case 6
                    D.trials(n).label = 'Fixation';
                otherwise
                    D.trials(n).label = 'Undefined';
            end
        end
    case 'StoryM'
    otherwise
end
D.path = savePath;
D.fname = ['spm_' fname '.mat'];

fnamedat = ['spm_' fname '.dat'];

if ~isTF
    if Ntrials == 1
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nsamples], 'float32-le');
        datafile(:, :) = data;
    else
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nsamples Ntrials], 'float32-le');
        datafile(:, :, :) = data;
    end
else
    if Ntrials == 1
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nfrequencies Nsamples], 'float32-le');
        datafile(:, :, :) = data;
    else
        datafile = file_array(fullfile(D.path, fnamedat), [Nchannels Nfrequencies Nsamples Ntrials], 'float32-le');
        datafile(:, :, :, :) = data;
    end    
    D.transform.ID = 'TF';
    D.transform.frequencies = ftdata.freq;
end

D.data = datafile;

D = meeg(D);

if  isfield(ftdata, 'hdr')
    % Uses fileio function to get the information about channel types stored in
    % the original header. This is now mainly useful for Neuromag support but might
    % have other functions in the future.
    origchantypes = ft_chantype(ftdata.hdr);
    [sel1, sel2] = spm_match_str(D.chanlabels, ftdata.hdr.label);
    origchantypes = origchantypes(sel2);
    if length(strmatch('unknown', origchantypes, 'exact')) ~= numel(origchantypes)
        D.origchantypes = struct([]);
        D.origchantypes(1).label = ftdata.hdr.label(sel2);
        D.origchantypes(1).type = origchantypes;
    end
end

% Set channel types to default
S1 = [];
S1.task = 'defaulttype';
S1.D = D;
S1.updatehistory = 0;
D = spm_eeg_prep(S1);

if Ntrials == 1
    D = type(D, 'continuous');
else
    D = type(D, 'single');
end

if  isfield(ftdata, 'hdr') && isfield(ftdata.hdr, 'grad')
    D = sensors(D, 'MEG', ft_convert_units(ftdata.hdr.grad, 'mm'));
    
    S = [];
    S.task = 'project3D';
    S.modality = 'MEG';
    S.updatehistory = 0;
    S.D = D;

    D = spm_eeg_prep(S);
end
presentDir = pwd;
cd(savePath);
D = check(D);
cd(presentDir);
save(D);
end


function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning a cell array
%
% Use as
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where str is a string and sep is the separator at which you want
% to cut it into pieces.
%
% Using the optional boolean flag rep you can specify whether repeated
% seperator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% seperators are treated as one.

% Copyright (C) 2003-2006, Robert Oostenveld


tok = {};
f = find(str==sep);
f = [0, f, length(str)+1];
for i=1:(length(f)-1)
    tok{i} = str((f(i)+1):(f(i+1)-1));
end

if nargin<3 || rep
    % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
    tok(find(cellfun('isempty', tok)))=[];
end
end
