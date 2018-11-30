%% Load Raw data
addpath('E:\my_code');
dataDir = 'E:\MouseLFP\11mo\WT';
workingDir = 'E:\LFP_Example\DCMs\Group';

dataDirInfo = dir(dataDir);
subject = [];
for i=1:length(dataDirInfo)
    if ~isempty(strfind(dataDirInfo(i).name, '#')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(dataDirInfo(i).name,'#');
        subject(end).dir = fullfile(dataDir,dataDirInfo(i).name);
    end
end
for i =1:length(subject)
    subjDirInfo = dir(subject(i).dir);
    subject(i).files = {};
    subject(i).ntrial = 0;
    for j = 1:length(subjDirInfo)
        if ~isempty(strtok(subjDirInfo(j).name,'.'))
            subject(i).files{end+1} = subjDirInfo(j).name;
            subject(i).ntrial = subject(i).ntrial +1;
        end
    end
end
clear subjDirInfo i j;

cd(workingDir);
dataset = cell(length(subject), 1);

%cfg.savePath                 = workingDir;

for i = 1:length(subject)
    for session = 1:length(subject(i).files)
        fid = fopen(fullfile(subject(i).dir,subject(i).files{session}));
        Raw_d{i,session} = textscan(fid, '%f %f  %f  %f %f  %f %f  %f',20000);
    end
    fclose('all');
end
clear i session
resultDIr = 'E:\LFP_Example\DCMs';
%% Check raw data and assign to fieldtrip data structures
for i = 1:length(subject)
    isdataok = 1;
    for session = 1:subject(i).ntrial
        if length(Raw_d{i,session})~=8
            isdataok = 0;
            error('This session has less than 8 channels')
        end
        for channel = 1:8
            if length(Raw_d{session}{channel})~=16000
                isdataok = 0;
                error ('This trial has more than 16000 time point');
            end
        end
    end
    if isdataok
        fprintf('Raw data of subject %s is successfully loaded.\n',subject(i).id);
       %% Create fieldtrip type data(Concatenate all trials)
        data = [];
        for session = 1:subject(i).ntrial
            start(session) = ((session-1)*16000+1);
            epoch = [start(session): (start(session)+16000-1)];
            for chan = 1:8
                data(epoch,chan) = Raw_d{i,session}{chan};
            end
        end
        
        ftdata = [];
        ftdata.trial{1} = data(:,1:8)';
        ftdata.time{1} = (0:(16000*subject(i).ntrial - 1))./1600;
        
        % Channel define
        % 1. right ACC  2.left ACC   3. right AI(Anterior Insula)   4. left AI
        % 5. right RSC(Retrosplenial cortex)  6. left RSC
        % 7. right BLA(Basolateral Amygdala   8. left BLA
        chanlabel ={'rACC','lACC','rAI', 'lAI', 'rRSC', 'lRSC', 'rBLA','lBLA'};
        ftdata.fsample = 1600;
        ftdata.label = chanlabel(1:8);
        ftdata.label = ftdata.label(:);
        dataset{i} = ftdata;
    end
end
clear start session i Raw_d isdataok fid ftdata epoch dataDirInfo data channel chan chanlabel
%% Preprocess(Downsample -> High pass filter-> Notch filter-> Low pass filter)
% Plot rawdata
dataclean = {};
for i = 1:length(subject)
    figure;
    subplot(6,2,1);
    plot(dataset{i}.time{1}, dataset{i}.trial{1});
    grid(gca,'on')
    xticks(gca, 10:10:10*subject(i).ntrial);
    xlim(gca, [0,10*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('raw data timeseries(concatenated all trials)');
    subplot(6,2,2);
    length(dataset{i}.time)
    plot((0:(length(dataset{i}.time{1})-1))*...
        (dataset{i}.fsample/length(dataset{i}.time{1})),abs(fft(dataset{i}.trial{1}(1,:))));
    title('raw data powerspectrum');
    
    % Downample to 800Hz
    cfg = [];
    cfg.resamplefs = 800;
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    datadown = ft_resampledata(cfg, dataset{i});
    
    subplot(6,2,3);
    plot(datadown.time{1}, datadown.trial{1});
    grid(gca,'on')
    xticks(gca, 10:10:10*subject(i).ntrial);
    xlim(gca, [0,10*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('downsampled data timeseries (800 Hz)');
    subplot(6,2,4);
    plot((0:(length(datadown.time{1})-1))*...
        (datadown.fsample/length(datadown.time{1})),abs(fft(datadown.trial{1}(1,:))));
    title('downsampled data powerspectrum');
    
    % Band pass filter (1, 200) Hz
    cfg = [];
    cfg.bpfilter ='yes';
    cfg.bpfreq = [0.5, 200]; %% must satisfy this condition ( freqmax*2 < fsample )
    databp = ft_preprocessing(cfg, datadown);
    
    subplot(6,2,5);
    plot(databp.time{1}, databp.trial{1});
    grid(gca,'on')
    xticks(gca, 10:10:10*subject(i).ntrial);
    xlim(gca, [0,10*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('Band pass filtered data timeseries (1-200 Hz)');
    subplot(6,2,6);
    plot((0:(length(databp.time{1})-1))*...
        (databp.fsample/length(databp.time{1})),abs(fft(databp.trial{1}(1,:))));
    title('Band pass filtered data powerspectrum');
    
    % Notch filter(60Hz, 120Hz, 180Hz)
    cfg = [];
    cfg.dftfreq   = [60-1:(1/100):60+1 120-1:(1/100):120+1 ...
        180-1:(1/100):180+1]; % filter out 60 hz line noise
    cfg.dftfilter = 'yes';
    datanotch     = ft_preprocessing(cfg,databp);
    
    subplot(6,2,7);
    plot(datanotch.time{1}, datanotch.trial{1});
    grid(gca,'on')
    xticks(gca, 10:10:10*subject(i).ntrial);
    xlim(gca, [0,10*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('Notch filtered data timeseries (60, 120, 180 Hz)');
    subplot(6,2,8);
    plot((0:(length(datanotch.time{1})-1))*...
        (datanotch.fsample/length(datanotch.time{1})),abs(fft(datanotch.trial{1}(1,:))));
    title('Notch pass filtered data powerspectrum');

    
    % Erase trial Edge effect by remove edges of each trials                              
    % for time window [trial_start, trial_start + 400ms], 
    % [trial_end -100 ms, trial_end].
    dataedge = datanotch;
    edgezoneidx = [];

    for j=0:(subject(i).ntrial)
        edgezoneidx = [edgezoneidx, find(dataedge.time{1}<j*10+0.5 & dataedge.time{1}>=j*10-0.5)];
    end
    dataedge.trial{1}(:,edgezoneidx) = [];
    dataedge.time{1} = [];
    dataedge.time{1} = ((1:length(dataedge.trial{1}))-1)*(1/dataedge.fsample);
    dataedge.sampleinfo = [1,length(dataedge.time{1})];
    
    subplot(6,2,9);
    plot(dataedge.time{1}, dataedge.trial{1});
    grid(gca,'on')
    xticks(gca, 9:9:9*subject(i).ntrial);
    xlim(gca, [0,9*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('Trial edge cleard data timeseries');
    subplot(6,2,10);
    plot((0:(length(dataedge.time{1})-1))*...
        (dataedge.fsample/length(dataedge.time{1})),abs(fft(dataedge.trial{1}(1,:))));
    title('Edge cleared data powerspectrum');
    
    
    % Calculate Mean and standard deviation, ceiling irregular spikes with
    % [mean - 3*stddev, mean + 3*stddev]
    dataclean{i} = dataedge;
    lowerbound = [];upperbound = [];
    lowerbound = (mean(dataclean{i}.trial{1}') - 3 * std(dataclean{i}.trial{1}'))';
    upperbound = (mean(dataclean{i}.trial{1}') + 3 * std(dataclean{i}.trial{1}'))';
    overpeak = dataclean{i}.trial{1} > upperbound;
    underpeak = dataclean{i}.trial{1} < lowerbound;
    dataclean{i}.trial{1} = dataclean{i}.trial{1}.*~overpeak + overpeak.*repmat(upperbound,1,length(dataclean{i}.time{1}));
    dataclean{i}.trial{1} = dataclean{i}.trial{1}.*~underpeak + underpeak.*repmat(lowerbound,1,length(dataclean{i}.time{1}));
    
    subplot(6,2,11);
    plot(dataclean{i}.time{1}, dataclean{i}.trial{1});
    grid(gca,'on')
    xticks(gca, 9:9:9*subject(i).ntrial);
    xlim(gca, [0,9*subject(i).ntrial]);
    ylim(gca, [-2000 2000]);
    title('Z score ceiling data timeseries(3 standard deviation)');
    subplot(6,2,12);
    plot((0:(length(dataclean{i}.time{1})-1))*...
        (dataclean{i}.fsample/length(dataclean{i}.time{1})),abs(fft(dataclean{i}.trial{1}(1,:))));
    title('Z score ceiling data powerspectrum');
    suptitle(['Subject ', subject(i).id]);
end
clear data cfg chanlabel databp datadown i j edgezoneidx lowerbound overpeak underpeak upperbound dataedge datanotch
%% Subject 84 has too much noisy trial(reject 18s - 27s).
dataidx = find(strcmp({subject.id},'84'));
figure;
subplot(2,1,1);
plot(dataclean{dataidx}.time{1}, dataclean{dataidx}.trial{1});
grid(gca,'on')
xticks(gca, 9:9:9*subject(dataidx).ntrial);
xlim(gca, [0,9*subject(dataidx).ntrial]);
ylim(gca, [-2000 2000]);
title('subject 84 before trial rejection');

mask = [];
mask = find(dataclean{dataidx}.time{1}>=18 & dataclean{dataidx}.time{1}<27);
dataclean{dataidx}.trial{1}(:,mask) = [];
dataclean{dataidx}.time{1} = [];
dataclean{dataidx}.time{1} = ((1:length(dataclean{dataidx}.trial{1}))-1)*(1/dataclean{dataidx}.fsample);
dataclean{dataidx}.sampleinfo = [1,length(dataclean{dataidx}.time{1})];

subplot(2,1,2);
plot(dataclean{dataidx}.time{1}, dataclean{dataidx}.trial{1});
grid(gca,'on')
xticks(gca, 9:9:9*(subject(dataidx).ntrial-1));
xlim(gca, [0,9*(subject(dataidx).ntrial-1)]);
ylim(gca, [-2000 2000]);
title('subject 84 after trial rejection');
suptitle(['Subject ', subject(dataidx).id]);

clear mask dataidx
%% Epoch Data
bilateral_dataclean = dataclean;
% Use right hemisphere
for i = 1:length(dataclean)
    newlabel = {};
    newlabel = dataclean{i}.label([2,4,6,8]);
    dataclean{i}.label = newlabel;
    dataclean{i}.trial{1}([1,3,4,7],:) = [];
end
clear newlabel
ftdata_concat = dataclean;
% separate each trials
for i = 1:length(dataclean)
    dataclean{i} = rmfield(dataclean{i}, 'sampleinfo');
    dataclean{i}.trial ={}; dataclean{i}.time = {};
    fsample = dataclean{i}.fsample;
    totlength = fsample * 1;
    time = (0:(totlength-1))./fsample;
    start = 1:totlength:length(ftdata_concat{i}.time{1});
    if ~~mod(length(ftdata_concat{i}.time{1}),totlength)
        ntrial = length(start) -1;
    else
        ntrial = length(start);
    end
    for trl = 1:ntrial
        dataind = (0:(totlength-1))+start(trl);
        dataclean{i}.trial{trl} = ftdata_concat{i}.trial{1}(:,dataind);
        dataclean{i}.time{trl} = time;
    end
end
ftdata = dataclean;
clear dataclean dataind fsample i ntrial start time totlength trl
%% Create spm MEEG datafiles
% define the output file name
for i = 1:length(ftdata)
    fname = ['left_LFP_WT_subj_' subject(i).id];
    % fname_concat = 'LFP_WT_subj4_concat';
    
    % Convert the ftdata struct to SPM M\EEG dataset
    %--------------------------------------------------------------------------
    D{i} = spm_eeg_ft2spm(ftdata{i}, fname);
    % D_concat = spm_eeg_ft2spm(ftdata_concat, fname_concat);
    
    % Examples of providing additional information in a script
    % ':' comes instead of an index vector and means that the command
    % applies to all channels/all trials.
    %--------------------------------------------------------------------------
    D{i} = type(D{i}, 'single');             % Sets the dataset type
    D{i} = chantype(D{i}, ':', 'LFP');        % Sets the channel type
    D{i} = conditions(D{i}, ':', 'WT-11');  % Sets the condition label
    
%     D_concat = type(D_concat, 'continuous');             % Sets the dataset type
%     D_concat = chantype(D_concat, ':', 'LFP');        % Sets the channel type
%     D_concat = conditions(D_concat, ':', 'WT-11');  % Sets the condition label
    % save
    %--------------------------------------------------------------------------
    save(D{i});
    %save(D_concat);
end
%% Run spDCM
if ~isdir('DCMs')
    mkdir('DCMs');
end

DCM = [];
DCM.xY.modality = 'LFP';
DCM.options.analysis = 'CSD';
DCM.options.model = 'CMC';
DCM.options.spatial = 'LFP';
DCM.options.trials = 1;
DCM.options.Tdcm = [1 1000]; % use whole time
DCM.options.Fdcm = [1 60];
DCM.options.Nmodes = 8;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
for i = 1:length(D)
    DCM.xY.Dfile = D{i}.fullfile;
    DCM.name = ['DCM_' D{i}.fname];
    %DCM.val =D.val;
    DCM_tmp = spm_dcm_erp_data(DCM, 0); % 0 for not to make erp
    DCM_tmp.Lpos   = [];
    DCM_tmp.Sname  = ftdata{i}.label;
    DCM_tmp.xU.X = [];
    DCM_tmp.xU.name = {};
    DCM_tmp.B ={};
    DCM_tmp.A{1} = ones(4);
    DCM_tmp.A{2} = ones(4);
    DCM_tmp.A{3} = ones(4);
    DCM_Result{i} = spm_dcm_csd(DCM_tmp);
end
DCM = DCM_Result;
clear DCM_tmp DCM_Result
%% PEB
M.X = ones(7,1);
PEBi = spm_dcm_peb(DCM,M);
for i = 1:7
    pebi{i,1} = PEBi(i);
end
PEBa = spm_dcm_peb_peb(pebi);