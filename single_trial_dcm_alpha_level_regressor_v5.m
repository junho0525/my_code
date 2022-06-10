%% get directories and subjects
addpath('E:\my_code'); 
clear;
cd E:\MEGResult\SingleTrialDCM\DCMs
dirInfo = dir();
isdirnames = zeros(1, length(dirInfo));
for i=1:length(dirInfo)
    if isfolder(dirInfo(i).name)&&~isnan(str2double(dirInfo(i).name(1)))
        isdirnames(i) = 1;
    end
end
dirInfo(isdirnames==0) = [];
clear isdirnames i
subjectnames = cell(1,length(dirInfo));
for i=1:length(dirInfo)
    subjectnames{i} = strtok(dirInfo(i).name,'_');
end
subjectnames = unique(subjectnames);
clear i dirInfo
%% 
use_my_library('spm',1);
mPath='E:\HCPMEG'; % This is the directory of root of HCP dataset.
workingDir = 'E:\MEGResult'; % This is the directory that all result files will be placed.
megDirInfo = dir(mPath);
subject = [];
cd(mPath);
for idx=1:length(megDirInfo)
    if ~isdir(megDirInfo(idx).name)
        continue;
    end
    if ~isempty(strfind(megDirInfo(idx).name, 'MEG')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(megDirInfo(idx).name,'_');
        subject(end).dir = megDirInfo(idx).name;
    elseif ~isempty(str2num(megDirInfo(idx).name)) % if subject folder is '[HCP data path]/[subject id]'
        subject(end+1).id = megDirInfo(idx).name;
        subject(end).dir = megDirInfo(idx).name;
    end
end
clear megDirInfo;
for idx =1:length(subject)
    subjDirInfo = dir(fullfile(mPath,subject(idx).dir));
    subject(idx).anatomy = [];
    subject(idx).Restin = [];
    subject(idx).Wrkmem = [];
    subject(idx).StoryM = [];
    subject(idx).Motort = [];
    for j =1:length(subjDirInfo)
       if ~isempty(strfind(subjDirInfo(j).name, 'anatomy'))
           subject(idx).anatomy = 1; 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Restin_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Restin/rmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Restin(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Restin =unique(subject(idx).Restin); 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Wrkmem_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Wrkmem/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Wrkmem(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Wrkmem =unique(subject(idx).Wrkmem);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'StoryM_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/StoryM/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).StoryM(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).StoryM =unique(subject(idx).StoryM);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Motort_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Motort/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Motort(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Motort =unique(subject(idx).Motort);
       end
    end
end
clear tempDirInfo subjDirInfo subDirInfo sessionNum i j k idx;
cd ('E:\MEGResult\SingleTrialDCM\DCMs');
isDCM = zeros(1,length(subject));
for idx = 1:length(subject)
    isDCM(idx) = sum(contains(subjectnames,subject(idx).id));
end
subject(~isDCM) = [];
clear isDCM idx

% random sampling of subjects( n = 5 )
subjectfilter = randperm(15);
subjectfilter(6:end) = [];
allsubject = subject;
subject = subject(subjectfilter);
allsubjectnames = subjectnames;
subjectnames = subjectnames(subjectfilter);
%% make connectivity models ( check every Full model - one connectivity )
A = {};
% full model
A{1}{1} = ones(4)*0.625;
A{1}{2} = ones(4)*0.625;
A{1}{3} = ones(4)*0.625;
for modelidx = 1:48
    A{modelidx+1} = A{1};
    A{modelidx+1}{floor((modelidx-1)/16)+1}((mod(modelidx-1,16)+1))= 0;
end
clear modelidx
%%
use_my_library('mnet',1);
designnames = {'Avg Pow', 'Alpha','Beta','Delta', 'Theta','Gamma'};
dcmcount = zeros(length(subject),1);
BMC2 = {};
pebm = {};
pebi = {};
for subjidx = 1:length(subject)
    %% set DCM names
    DCMs = {};
    DCMstemp = {};
    DCMlength = zeros(3,1);
    for sessionidx = 1:3
        subjectInfo = dir([subjectnames{subjidx} '_' num2str(sessionidx+2)]);
        isDCM = zeros(1,length(subjectInfo));
        for idx=1:length(subjectInfo)
            isDCM(idx) = contains(subjectInfo(idx).name, 'DCM');
        end
        subjectInfo(~isDCM) = [];
        clear idx isDCM
        for idx = 1:length(subjectInfo)
            if isfile(fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']))
                DCMs{idx,sessionidx} = fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']);
            else
                error(sprintf('There is no DCM file! : %s',fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat'])));
            end
        end
        DCMlength(sessionidx) =length(subjectInfo);
    end
    clear subjectInfo sessionidx idx
    %% Select DCMs according to Free Energy (Positive Free Enegy
    dcmFreeEnergy = {};
    dcmFreeEnergy{1}= zeros(DCMlength(1),1);
    dcmFreeEnergy{2}= zeros(DCMlength(2),1);
    dcmFreeEnergy{3}= zeros(DCMlength(3),1);
    for sessionidx =1:3
        for idx = 1:DCMlength(sessionidx)
            DCM = load(DCMs{idx,sessionidx});
            DCM = DCM.DCM;
            dcmFreeEnergy{sessionidx}(idx) = DCM.F;
        end
    end
    fprintf('Data Load End\n');
    dcmfilter = {};
    dcmfilter{1} = find(dcmFreeEnergy{1}>0);
    dcmfilter{2} = find(dcmFreeEnergy{2}>0);
    dcmfilter{3} = find(dcmFreeEnergy{3}>0);
    dcmcount(subjidx) = length(dcmfilter{1})+length(dcmfilter{2})+length(dcmfilter{3});
    DCMstemp(1:length(dcmfilter{1}),1) = DCMs(dcmfilter{1},1);
    DCMstemp(1:length(dcmfilter{2}),2) = DCMs(dcmfilter{2},2);
    DCMstemp(1:length(dcmfilter{3}),3) = DCMs(dcmfilter{3},3);
    DCMs = DCMstemp;
    clear DCMstemp dcmFreeEnergy dcmcount idx;
    DCMlength = [length(dcmfilter{1}), length(dcmfilter{2}),length(dcmfilter{3})];
    
    %% Compute band power and create design matrix
    alphalevels = {};
    averagepowlevels = {};
    betalevels = {};
    deltalevels = {};
    thetalevels = {};
    gammalevels = {};
    %figure;
    for sessionidx = 1:3
        data = [];
        datapow  = [];
        load (fullfile(mPath, subject(subjidx).dir, [subject(subjidx).dir '_Restin_preproc'], subject(subjidx).id, 'MEG/Restin/rmegpreproc' , ...
            [subject(subjidx).dir '_' num2str(sessionidx+2) '-Restin_rmegpreproc.mat']));
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper ='dpss';
        cfg.tapsmofrq = 1;
        cfg.keeptrials = 'yes';
        
        cfg.foilim = [8,12];
        datapow = ft_freqanalysis(cfg, data);
        alphalevels{sessionidx} = sum(sum(datapow.powspctrm,3),2);
        alphalevels{sessionidx} = alphalevels{sessionidx}(dcmfilter{sessionidx});
    end
    clear data datapow dcmfilter cfg
    
    designmatrix = {};
    for sessionidx = 1:3
        designmatrix{sessionidx}(:,1) = alphalevels{sessionidx};
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) - mean(designmatrix{sessionidx}(:,:));
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) ./ sqrt(var(designmatrix{sessionidx}(:,:)));
    end
    M.X = [];
    M.X = ones(sum(DCMlength),1);
    
    M.X(1:DCMlength(1),2) = designmatrix{1};
    M.X(DCMlength(1)+1:DCMlength(1)+DCMlength(2),2) = designmatrix{2};
    M.X(DCMlength(1)+DCMlength(2)+1:sum(DCMlength),2) = designmatrix{3};
    M.X(:,2:end) = M.X(:,2:end) - mean(M.X(:,2:end));
    M.X(:,2:end) = M.X(:,2:end)./sqrt(var(M.X(:,2:end)));
    
    clear alphalevels averagepowlevels betalevels thetalevels gammalevels deltalevels sessionidx designmatrix
    %% Specify Reduced model
    DCM = load(DCMs{1,1});
    DCM = DCM.DCM;
    RCM = {};
    RCM{1}.M.pE = DCM.M.pE; % RFCM{1} : fullmodel
    RCM{1}.M.pC = DCM.M.pC;
    
    for i=1:length(A) % number of reduced models
        RCM{i}.M.pE = DCM.M.pE;
        RCM{i}.M.pC = DCM.M.pC;
        RCM{i}.M.pC.A = A{i};
    end
    
    DCMall = {};
    DCMall = DCMs(1:DCMlength(1),1);
    DCMall(DCMlength(1)+1:DCMlength(1)+DCMlength(2)) = DCMs(1:DCMlength(2),2);
    DCMall(DCMlength(1)+DCMlength(2)+1:sum(DCMlength)) = DCMs(1:DCMlength(3),3);
    
    dcmlength = length(DCMall);
    
    for i = 1:dcmlength
        DCMall(i,2:length(A)) = RCM(2:length(A));
    end
    %% Bayesian Model Reduction -> Parametric Empirical Bayes(BMC_PEB)
    rcm = {};
    rcmtemp = [];
    BMC = [];
    % BMR
    [rcmtemp, BMC,~] = spm_dcm_bmr(DCMall,'A');
    for i = 1:length(rcmtemp)
        rcm(i,1:length(A)) = rcmtemp{i};
    end
    clear rcmtemp;
    % BMC_PEB
    for modelidx = 1:size(rcm,2)
        [BMC2{subjidx,modelidx},pebm{subjidx,modelidx}] = spm_dcm_bmc_peb(rcm(:,modelidx),M);
        pebi{subjidx, modelidx} = spm_dcm_peb(rcm(:,modelidx),M);
        drawnow;
    end
end
F = [];
peba = [];
for modelidx = 1:size(pebi,2)
    peba{modelidx} = mnet_dcm_peb_peb(pebi(:,modelidx));
    F(modelidx) = peba{modelidx}.F;
end
figure;
bar(F);
fprintf('Finish!\n');
%% make connectivity models ( check every refmodel + one connectivity )
A = {};
% full model
A{1}{1} = [0,0,0,0;0,0,0,0;1,1,0,0;1,1,1,0]*0.625;
A{1}{2} = A{1}{1}';
A{1}{3} = [0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]*0.625;
for modelidx = 1:48
    A{modelidx+1} = A{1};
    A{modelidx+1}{floor((modelidx-1)/16)+1}((mod(modelidx-1,16)+1))= 0.625;
end
clear modelidx
%%
use_my_library('mnet',1);
designnames = {'Avg Pow', 'Alpha','Beta','Delta', 'Theta','Gamma'};
dcmcount = zeros(length(subject),1);
BMC2 = {};
pebm = {};
pebi = {};
for subjidx = 1:length(subject)
    %% set DCM names
    DCMs = {};
    DCMstemp = {};
    DCMlength = zeros(3,1);
    for sessionidx = 1:3
        subjectInfo = dir([subjectnames{subjidx} '_' num2str(sessionidx+2)]);
        isDCM = zeros(1,length(subjectInfo));
        for idx=1:length(subjectInfo)
            isDCM(idx) = contains(subjectInfo(idx).name, 'DCM');
        end
        subjectInfo(~isDCM) = [];
        clear idx isDCM
        for idx = 1:length(subjectInfo)
            if isfile(fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']))
                DCMs{idx,sessionidx} = fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']);
            else
                error(sprintf('There is no DCM file! : %s',fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat'])));
            end
        end
        DCMlength(sessionidx) =length(subjectInfo);
    end
    clear subjectInfo sessionidx idx
    %% Select DCMs according to Free Energy (Positive Free Enegy
    dcmFreeEnergy = {};
    dcmFreeEnergy{1}= zeros(DCMlength(1),1);
    dcmFreeEnergy{2}= zeros(DCMlength(2),1);
    dcmFreeEnergy{3}= zeros(DCMlength(3),1);
    for sessionidx =1:3
        for idx = 1:DCMlength(sessionidx)
            DCM = load(DCMs{idx,sessionidx});
            DCM = DCM.DCM;
            dcmFreeEnergy{sessionidx}(idx) = DCM.F;
        end
    end
    fprintf('Data Load End\n');
    dcmfilter = {};
    dcmfilter{1} = find(dcmFreeEnergy{1}>0);
    dcmfilter{2} = find(dcmFreeEnergy{2}>0);
    dcmfilter{3} = find(dcmFreeEnergy{3}>0);
    dcmcount(subjidx) = length(dcmfilter{1})+length(dcmfilter{2})+length(dcmfilter{3});
    DCMstemp(1:length(dcmfilter{1}),1) = DCMs(dcmfilter{1},1);
    DCMstemp(1:length(dcmfilter{2}),2) = DCMs(dcmfilter{2},2);
    DCMstemp(1:length(dcmfilter{3}),3) = DCMs(dcmfilter{3},3);
    DCMs = DCMstemp;
    clear DCMstemp dcmFreeEnergy dcmcount idx;
    DCMlength = [length(dcmfilter{1}), length(dcmfilter{2}),length(dcmfilter{3})];
    
    %% Compute band power and create design matrix
    alphalevels = {};
    averagepowlevels = {};
    betalevels = {};
    deltalevels = {};
    thetalevels = {};
    gammalevels = {};
    %figure;
    for sessionidx = 1:3
        data = [];
        datapow  = [];
        load (fullfile(mPath, subject(subjidx).dir, [subject(subjidx).dir '_Restin_preproc'], subject(subjidx).id, 'MEG/Restin/rmegpreproc' , ...
            [subject(subjidx).dir '_' num2str(sessionidx+2) '-Restin_rmegpreproc.mat']));
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper ='dpss';
        cfg.tapsmofrq = 1;
        cfg.keeptrials = 'yes';
        
        cfg.foilim = [8,12];
        datapow = ft_freqanalysis(cfg, data);
        alphalevels{sessionidx} = sum(sum(datapow.powspctrm,3),2);
        alphalevels{sessionidx} = alphalevels{sessionidx}(dcmfilter{sessionidx});
    end
    clear data datapow dcmfilter cfg
    
    designmatrix = {};
    for sessionidx = 1:3
        designmatrix{sessionidx}(:,1) = alphalevels{sessionidx};
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) - mean(designmatrix{sessionidx}(:,:));
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) ./ sqrt(var(designmatrix{sessionidx}(:,:)));
    end
    M.X = [];
    M.X = ones(sum(DCMlength),1);
    
    M.X(1:DCMlength(1),2) = designmatrix{1};
    M.X(DCMlength(1)+1:DCMlength(1)+DCMlength(2),2) = designmatrix{2};
    M.X(DCMlength(1)+DCMlength(2)+1:sum(DCMlength),2) = designmatrix{3};
    M.X(:,2:end) = M.X(:,2:end) - mean(M.X(:,2:end));
    M.X(:,2:end) = M.X(:,2:end)./sqrt(var(M.X(:,2:end)));
    
    clear alphalevels averagepowlevels betalevels thetalevels gammalevels deltalevels sessionidx designmatrix
    %% Specify Reduced model
    DCM = load(DCMs{1,1});
    DCM = DCM.DCM;
    RCM = {};
    RCM{1}.M.pE = DCM.M.pE; % RFCM{1} : fullmodel
    RCM{1}.M.pC = DCM.M.pC;
    
    for i=1:length(A) % number of reduced models
        RCM{i}.M.pE = DCM.M.pE;
        RCM{i}.M.pC = DCM.M.pC;
        RCM{i}.M.pC.A = A{i};
    end
    
    DCMall = {};
    DCMall = DCMs(1:DCMlength(1),1);
    DCMall(DCMlength(1)+1:DCMlength(1)+DCMlength(2)) = DCMs(1:DCMlength(2),2);
    DCMall(DCMlength(1)+DCMlength(2)+1:sum(DCMlength)) = DCMs(1:DCMlength(3),3);
    
    dcmlength = length(DCMall);
    
    for i = 1:dcmlength
        DCMall(i,2:length(A)) = RCM(2:length(A));
    end
    %% Bayesian Model Reduction -> Parametric Empirical Bayes(BMC_PEB)
    rcm = {};
    rcmtemp = [];
    BMC = [];
    % BMR
    [rcmtemp, BMC,~] = spm_dcm_bmr(DCMall,'A');
    for i = 1:length(rcmtemp)
        rcm(i,1:length(A)) = rcmtemp{i};
    end
    clear rcmtemp;
    % BMC_PEB
    for modelidx = 1:size(rcm,2)
        [BMC2{subjidx,modelidx},pebm{subjidx,modelidx}] = spm_dcm_bmc_peb(rcm(:,modelidx),M);
        pebi{subjidx, modelidx} = spm_dcm_peb(rcm(:,modelidx),M);
        drawnow;
    end
end
F2 = [];
peba = [];
for modelidx = 1:size(pebi,2)
    peba{modelidx} = mnet_dcm_peb_peb(pebi(:,modelidx));
    F2(modelidx) = peba{modelidx}.F;
end
figure;
bar(F2);
fprintf('Finish!\n');
%% make connectivity models ( check reduced models of ref model )
A = {};
% ref model
A{1}{1} = [0,0,0,0;0,0,0,0;1,1,0,0;1,1,1,0]*0.625;
A{1}{2} = A{1}{1}';
A{1}{3} = [0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]*0.625;
A{2} = A{1};
clear modelidx
%%
use_my_library('mnet',1);
designnames = {'Avg Pow', 'Alpha','Beta','Delta', 'Theta','Gamma'};
dcmcount = zeros(length(subject),1);
BMC2 = {};
pebm = {};
pebi = {};
for subjidx = 1:length(subject)
    %% set DCM names
    DCMs = {};
    DCMstemp = {};
    DCMlength = zeros(3,1);
    for sessionidx = 1:3
        subjectInfo = dir([subjectnames{subjidx} '_' num2str(sessionidx+2)]);
        isDCM = zeros(1,length(subjectInfo));
        for idx=1:length(subjectInfo)
            isDCM(idx) = contains(subjectInfo(idx).name, 'DCM');
        end
        subjectInfo(~isDCM) = [];
        clear idx isDCM
        for idx = 1:length(subjectInfo)
            if isfile(fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']))
                DCMs{idx,sessionidx} = fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']);
            else
                error(sprintf('There is no DCM file! : %s',fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat'])));
            end
        end
        DCMlength(sessionidx) =length(subjectInfo);
    end
    clear subjectInfo sessionidx idx
    %% Select DCMs according to Free Energy (Positive Free Enegy
    dcmFreeEnergy = {};
    dcmFreeEnergy{1}= zeros(DCMlength(1),1);
    dcmFreeEnergy{2}= zeros(DCMlength(2),1);
    dcmFreeEnergy{3}= zeros(DCMlength(3),1);
    for sessionidx =1:3
        for idx = 1:DCMlength(sessionidx)
            DCM = load(DCMs{idx,sessionidx});
            DCM = DCM.DCM;
            dcmFreeEnergy{sessionidx}(idx) = DCM.F;
        end
    end
    fprintf('Data Load End\n');
    dcmfilter = {};
    dcmfilter{1} = find(dcmFreeEnergy{1}>0);
    dcmfilter{2} = find(dcmFreeEnergy{2}>0);
    dcmfilter{3} = find(dcmFreeEnergy{3}>0);
    dcmcount(subjidx) = length(dcmfilter{1})+length(dcmfilter{2})+length(dcmfilter{3});
    DCMstemp(1:length(dcmfilter{1}),1) = DCMs(dcmfilter{1},1);
    DCMstemp(1:length(dcmfilter{2}),2) = DCMs(dcmfilter{2},2);
    DCMstemp(1:length(dcmfilter{3}),3) = DCMs(dcmfilter{3},3);
    DCMs = DCMstemp;
    clear DCMstemp dcmFreeEnergy dcmcount idx;
    DCMlength = [length(dcmfilter{1}), length(dcmfilter{2}),length(dcmfilter{3})];
    
    %% Compute band power and create design matrix
    alphalevels = {};
    %figure;
    for sessionidx = 1:3
        data = [];
        datapow  = [];
        load (fullfile(mPath, subject(subjidx).dir, [subject(subjidx).dir '_Restin_preproc'], subject(subjidx).id, 'MEG/Restin/rmegpreproc' , ...
            [subject(subjidx).dir '_' num2str(sessionidx+2) '-Restin_rmegpreproc.mat']));
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper ='dpss';
        cfg.tapsmofrq = 1;
        cfg.keeptrials = 'yes';
        
        cfg.foilim = [8,12];
        datapow = ft_freqanalysis(cfg, data);
        alphalevels{sessionidx} = sum(sum(datapow.powspctrm,3),2);
        alphalevels{sessionidx} = alphalevels{sessionidx}(dcmfilter{sessionidx});
    end
    clear data datapow dcmfilter cfg
    
    designmatrix = {};
    for sessionidx = 1:3
        designmatrix{sessionidx}(:,1) = alphalevels{sessionidx};
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) - mean(designmatrix{sessionidx}(:,:));
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) ./ sqrt(var(designmatrix{sessionidx}(:,:)));
    end
    M.X = [];
    M.X = ones(sum(DCMlength),1);
    
    M.X(1:DCMlength(1),2) = designmatrix{1};
    M.X(DCMlength(1)+1:DCMlength(1)+DCMlength(2),2) = designmatrix{2};
    M.X(DCMlength(1)+DCMlength(2)+1:sum(DCMlength),2) = designmatrix{3};
    M.X(:,2:end) = M.X(:,2:end) - mean(M.X(:,2:end));
    M.X(:,2:end) = M.X(:,2:end)./sqrt(var(M.X(:,2:end)));
    
    clear alphalevels averagepowlevels betalevels thetalevels gammalevels deltalevels sessionidx designmatrix
    %% Specify Reduced model
    DCM = load(DCMs{1,1});
    DCM = DCM.DCM;
    RCM = {};
    RCM{1}.M.pE = DCM.M.pE; % RFCM{1} : fullmodel
    RCM{1}.M.pC = DCM.M.pC;
    
    for i=1:length(A) % number of reduced models
        RCM{i}.M.pE = DCM.M.pE;
        RCM{i}.M.pC = DCM.M.pC;
        RCM{i}.M.pC.A = A{i};
    end
    
    DCMall = {};
    DCMall = DCMs(1:DCMlength(1),1);
    DCMall(DCMlength(1)+1:DCMlength(1)+DCMlength(2)) = DCMs(1:DCMlength(2),2);
    DCMall(DCMlength(1)+DCMlength(2)+1:sum(DCMlength)) = DCMs(1:DCMlength(3),3);
    
    dcmlength = length(DCMall);
    
    for i = 1:dcmlength
        DCMall(i,2:length(A)) = RCM(2:length(A));
    end
    %% Bayesian Model Reduction -> Parametric Empirical Bayes(BMC_PEB)
    rcm = {};
    rcmtemp = [];
    BMC = [];
    % BMR
    [rcmtemp, BMC,~] = spm_dcm_bmr(DCMall,'A');
    for i = 1:length(rcmtemp)
        rcm(i,1:length(A)) = rcmtemp{i};
    end
    clear rcmtemp;
    % BMC_PEB
    for modelidx = 1:size(rcm,2)
        [BMC2{subjidx,modelidx},pebm{subjidx,modelidx}] = spm_dcm_bmc_peb(rcm(:,modelidx),M);
        pebi{subjidx, modelidx} = spm_dcm_peb(rcm(:,modelidx),M);
        drawnow;
    end
end
F2 = [];
peba = [];
for modelidx = 1:size(pebi,2)
    peba{modelidx} = mnet_dcm_peb_peb(pebi(:,modelidx));
    F2(modelidx) = peba{modelidx}.F;
end
figure;
bar(F2);
fprintf('Finish!\n');
%% 20190329 make connectivity models ( check ref models without any regressor )
A = {};
% ref model
% mPFC > Prec > rLP = lLP
A{1}{1} = [0,0,0,0;0,0,0,0;1,1,0,0;1,1,1,0]*0.625;
A{1}{2} = A{1}{1}';
A{1}{3} = [0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]*0.625;

% Prec > mPFC > rLP = lLP
A{2}{1} = [0,0,0,0;0,0,0,0;1,1,0,1;1,1,0,0]*0.625;
A{2}{2} = A{2}{1}';
A{2}{3} = [0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]*0.625;

% mPFC = Prec > rLP = lLP
A{3}{1} = [0,0,0,0;0,0,0,0;1,1,0,0;1,1,0,0]*0.625;
A{3}{2} = A{3}{1}';
A{3}{3} = [0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0]*0.625;

% mPFC > Prec = rLP = lLP
A{4}{1} = [0,0,0,0;0,0,0,0;0,0,0,0;1,1,1,0]*0.625;
A{4}{2} = A{4}{1}';
A{4}{3} = [0,1,1,0;1,0,1,0;1,1,0,0;0,0,0,0]*0.625;

% rLP = lLP > Prec > mPFC
A{5}{1} = [0,0,1,1;0,0,1,1;0,0,0,1;0,0,0,0]*0.625;
A{5}{2} = A{5}{1}';
A{5}{3} = [0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]*0.625;

% mPFC = Prec = lLP =rLP
A{6}{1} = zeros(4);
A{6}{2} = A{6}{1}';
A{6}{3} = ones(4)*0.625;

%%
use_my_library('mnet',1);
designnames = {'Avg Pow', 'Alpha','Beta','Delta', 'Theta','Gamma'};
dcmcount = zeros(length(subject),1);
BMC2 = {};
pebm = {};
pebi = {};
for subjidx = 1:length(subject)
    %% set DCM names
    DCMs = {};
    DCMstemp = {};
    DCMlength = zeros(3,1);
    for sessionidx = 1:3
        subjectInfo = dir([subjectnames{subjidx} '_' num2str(sessionidx+2)]);
        isDCM = zeros(1,length(subjectInfo));
        for idx=1:length(subjectInfo)
            isDCM(idx) = contains(subjectInfo(idx).name, 'DCM');
        end
        subjectInfo(~isDCM) = [];
        clear idx isDCM
        for idx = 1:length(subjectInfo)
            if isfile(fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']))
                DCMs{idx,sessionidx} = fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat']);
            else
                error(sprintf('There is no DCM file! : %s',fullfile([subjectnames{subjidx} '_' num2str(sessionidx+2)],...
                ['DCM_' subjectnames{subjidx} '_sess' num2str(sessionidx+2) '_trial' num2str(idx) '.mat'])));
            end
        end
        DCMlength(sessionidx) =length(subjectInfo);
    end
    clear subjectInfo sessionidx idx
    %% Select DCMs according to Free Energy (Positive Free Enegy
    dcmFreeEnergy = {};
    dcmFreeEnergy{1}= zeros(DCMlength(1),1);
    dcmFreeEnergy{2}= zeros(DCMlength(2),1);
    dcmFreeEnergy{3}= zeros(DCMlength(3),1);
    for sessionidx =1:3
        for idx = 1:DCMlength(sessionidx)
            DCM = load(DCMs{idx,sessionidx});
            DCM = DCM.DCM;
            dcmFreeEnergy{sessionidx}(idx) = DCM.F;
        end
    end
    fprintf('Data Load End\n');
    dcmfilter = {};
    dcmfilter{1} = find(dcmFreeEnergy{1}>0);
    dcmfilter{2} = find(dcmFreeEnergy{2}>0);
    dcmfilter{3} = find(dcmFreeEnergy{3}>0);
    dcmcount(subjidx) = length(dcmfilter{1})+length(dcmfilter{2})+length(dcmfilter{3});
    DCMstemp(1:length(dcmfilter{1}),1) = DCMs(dcmfilter{1},1);
    DCMstemp(1:length(dcmfilter{2}),2) = DCMs(dcmfilter{2},2);
    DCMstemp(1:length(dcmfilter{3}),3) = DCMs(dcmfilter{3},3);
    DCMs = DCMstemp;
    clear DCMstemp dcmFreeEnergy dcmcount idx;
    DCMlength = [length(dcmfilter{1}), length(dcmfilter{2}),length(dcmfilter{3})];
    
    %% Compute band power and create design matrix
    alphalevels = {};
    %figure;
    for sessionidx = 1:3
        data = [];
        datapow  = [];
        load (fullfile(mPath, subject(subjidx).dir, [subject(subjidx).dir '_Restin_preproc'], subject(subjidx).id, 'MEG/Restin/rmegpreproc' , ...
            [subject(subjidx).dir '_' num2str(sessionidx+2) '-Restin_rmegpreproc.mat']));
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper ='dpss';
        cfg.tapsmofrq = 1;
        cfg.keeptrials = 'yes';
        
        cfg.foilim = [8,12];
        datapow = ft_freqanalysis(cfg, data);
        alphalevels{sessionidx} = sum(sum(datapow.powspctrm,3),2);
        alphalevels{sessionidx} = alphalevels{sessionidx}(dcmfilter{sessionidx});
    end
    clear data datapow cfg dcmfilter
    
    designmatrix = {};
    for sessionidx = 1:3
        designmatrix{sessionidx}(:,1) = alphalevels{sessionidx};
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) - mean(designmatrix{sessionidx}(:,:));
        designmatrix{sessionidx}(:,:) = designmatrix{sessionidx}(:,:) ./ sqrt(var(designmatrix{sessionidx}(:,:)));
    end
    M.X = [];
    M.X = ones(sum(DCMlength),1);
    
    M.X(1:DCMlength(1),2) = designmatrix{1};
    M.X(DCMlength(1)+1:DCMlength(1)+DCMlength(2),2) = designmatrix{2};
    M.X(DCMlength(1)+DCMlength(2)+1:sum(DCMlength),2) = designmatrix{3};
    M.X(:,2:end) = M.X(:,2:end) - mean(M.X(:,2:end));
    M.X(:,2:end) = M.X(:,2:end)./sqrt(var(M.X(:,2:end)));
    
    
    clear alphalevels averagepowlevels betalevels thetalevels gammalevels deltalevels sessionidx designmatrix
    %% Specify Reduced model
    DCM = load(DCMs{1,1});
    DCM = DCM.DCM;
    RCM = {};
    RCM{1}.M.pE = DCM.M.pE; % RFCM{1} : fullmodel
    RCM{1}.M.pC = DCM.M.pC;
    
    for i=1:length(A) % number of reduced models
        RCM{i}.M.pE = DCM.M.pE;
        RCM{i}.M.pE.A{1} = ~A{i}{1}*-4;
        RCM{i}.M.pE.A{2} = ~A{i}{2}*-4;
        RCM{i}.M.pE.A{3} = ~A{i}{3}*-4;
        RCM{i}.M.pC = DCM.M.pC;
        RCM{i}.M.pC.A = A{i};
    end
    
    % First model is full model
    DCMall = {};
    DCMall = DCMs(1:DCMlength(1),1);
    DCMall(DCMlength(1)+1:DCMlength(1)+DCMlength(2)) = DCMs(1:DCMlength(2),2);
    DCMall(DCMlength(1)+DCMlength(2)+1:sum(DCMlength)) = DCMs(1:DCMlength(3),3);
    
    dcmlength = length(DCMall);
    
    for i = 1:dcmlength
        DCMall(i,(1:length(A))+1) = RCM(1:length(A));
    end
    %% Bayesian Model Reduction -> Parametric Empirical Bayes(BMC_PEB)
    rcm = {};
    rcmtemp = [];
    BMC = [];
    % BMR
    [rcmtemp, BMC,~] = spm_dcm_bmr(DCMall,'A');
%     for i = 1:length(rcmtemp)
%         rcm(i,1:length(rcmtemp{i})) = rcmtemp(:,i);
%     end
    rcm = rcmtemp;
    clear rcmtemp;
    % BMC_PEB
    for modelidx = 1:size(rcm,2)
        pebi{subjidx, modelidx} = spm_dcm_peb(rcm(:,modelidx),M);
        drawnow;
    end
end
F2 = [];
peba = [];
for modelidx = 1:size(pebi,2)
    peba{modelidx} = mnet_dcm_peb_peb(pebi(:,modelidx));
    F2(modelidx) = peba{modelidx}.F;
end
figure;
bar(F2);
fprintf('Finish!\n');


%%
F = [];

peba = [];
M2.X = ones(10,1);
peba = mnet_dcm_peb_peb(pebm,M2);

spm_plot_ci2(pebi{designidx,1}.Ep(:,1),pebi{designidx,1}.Cp(1:parameternum,1:parameternum));

for i = 1:7
    figure;
    for j = 1:3
        cind = (1:16) + (i-1)*48+(j-1)*16;
        subplot(1,3,j);
        mnet_dcm_adjacency_ci(full(reshape(peba.Ep((1:16)+16*(j-1),i),4,4)),full(peba.Cp(cind,cind)));
    end
    if i ==1
        suptitle('Average');
    else
        suptitle(designnames{i-1});
    end
end
pebb{1} = mnet_dcm_peb_peb(pebi,M2);
pebb{2} = mnet_dcm_peb_peb(pebnoreg,M2);
pebb{3} = mnet_dcm_peb_peb(pebavg,M2);
pebb{4} = mnet_dcm_peb_peb(pebnoavg,M2);

pebb{3} = mnet_dcm_peb_peb(pebalpha,M2);
pebb{4} = mnet_dcm_peb_peb(pebbeta,M2);
pebb{5} = mnet_dcm_peb_peb(pebdelta,M2);
pebb{6} = mnet_dcm_peb_peb(pebtheta,M2);
pebb{7} = mnet_dcm_peb_peb(pebgamma,M2);
F(1) = peba.F;

for i = 1:7
    F(i+1) = pebb{i}.F
end
figure;
bar(F);
title('bmc\_peb vs manual selection of regressor');
xticklabels({'bmc\_peb(full)','full reg','no regressor','alpha-', 'beta-','delta-','theta-','gamma-'});
for i = 1:7
    figure('Position',[200,200,1500,500]);
    for j = 1:3
        cind = (1:16) + (i-1)*48+(j-1)*16;
        subplot(1,3,j);
        mnet_dcm_adjacency_ci(full(reshape(pebb.Ep((1:16)+16*(j-1),i),4,4)),full(pebb.Cp(cind,cind)));
    end
    if i ==1
        suptitle('Average');
    else
        suptitle(designnames{i-1});
    end
end

