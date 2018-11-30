%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP MEG Source Level FC Pipeline Using Field Trip                       %
%                                                                         %
% This pipeline code sould be executed in -nodeskop mode                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.04.04 15:55 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fieldtrip Resting Source Level FC Analysis
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg = [];
cfg.taskType = 'Restin'; % Default
cfg.savePath = restWorkingDir;
cfg.showMEGSignal            = 0;
cfg.showForwardResult        = 0; % yes : 1, no : 0(default)
cfg.sourceLocalization       = 1;
cfg.showSensPower            = 0;
cfg.showAlphaPowRatio        = 0;
cfg.showConnectionDensityMap = 0;
cfg.showFullFC               = 0;
cfg.showSeedBasedFC          = 1; % yes : 1(default) no : 0
cfg.saveFig                  = 1;
for i = 1:length(subject)
    cfg.subjectID           = subject(i).id;
    cfg.subjectPath         = fullfile(mPath, subject(i).dir);
    cfg.startIdx            = subject(i).startIdx;
    for sess = subject(i).Restin
        cfg.sessionID       = sess;
        D = rMEG_pipeline_HCP_jhs_v5(cfg);
        cfg.startIdx        = 1;
        close all;
    end
end