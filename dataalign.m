
savePath = restWorkingDir;
rawData = [];
preprocessData = [];
for idx = 1:length(subject)
    spmSession = cell(length(subject(idx).Restin),1);
    fieldSession = cell(length(subject(idx).Restin),1);;
    for sess = 1:length(subject(idx).Restin)
        subjectID = subject(idx).id; 
        subjectPath = fullfile(mPath, subject(idx).dir);
        sessionID = subject(idx).Restin(sess);
        rMEGPath = fullfile(subjectPath,[subjectID '_MEG_Restin_preproc'],subjectID,'MEG','Restin','rmegpreproc');
        rMEGRawData=[subjectID '_MEG_' num2str(sessionID) '-Restin_rmegpreproc'];
        try 
            D = spm_eeg_load(fullfile(savePath,['affdspm_', rMEGRawData]));
            fieldtripData = load(fullfile(rMEGPath, rMEGRawData));
        end
        spmSession{sess} = D(:,:,:);
        tmp = [];
        for trial = 1:length(fieldtripData.data.trial)
            tmp(:,:,trial) = fieldtripData.data.trial{trial};
        end
        fieldSession{sess} = tmp;
    end
    if ~isempty(spmSession{sess})
        preprocessData(end+1).subjectID = subjectID;
        preprocessData(end).session =  spmSession;
        preprocessData(end).samplingRate = D.fsample;
    end
    if ~isempty(fieldSession{sess})
        rawData(end+1).subjectID = subjectID;
        rawData(end).session = fieldSession;
        rawData(end).samplingRate = fieldtripData.data.fsample;
    end
end