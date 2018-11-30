
savePath = restWorkingDir;

dataset = [];
for idx = 1:length(subject)
    spmSession = cell(length(subject(idx).Restin),1);
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
        
    end
    if ~isempty(spmSession{sess})
        dataset(end+1) = struct('subjectID', subjectID, 'preprecessData', spmSession,'rawData');
    end
end