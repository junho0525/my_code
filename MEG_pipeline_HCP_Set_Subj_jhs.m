%% Get Subject from HCP Dataset
use_my_library('ft',1);
clear;
mPath='/media/kuro/DATA1/HCPMEG';
megDirInfo = dir(mPath);
subject = [];
for i=1:length(megDirInfo)
    if ~isempty(strfind(megDirInfo(i).name, 'MEG')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(megDirInfo(i).name,'_');
        subject(end).dir = megDirInfo(i).name;
    elseif ~isempty(str2num(megDirInfo(i).name)) % if subject folder is '[HCP data path]/[subject id]'
        subject(end+1).id = megDirInfo(i).name;
        subject(end).dir = megDirInfo(i).name;
    end
end
clear megDirInfo;
for i =1:length(subject)
    subjDirInfo = dir(fullfile(mPath,subject(i).dir));
    subject(i).anatomy = [];
    subject(i).Restin = [];
    subject(i).Wrkmem = [];
    subject(i).StoryM = [];
    subject(i).Motort = [];
    subject(i).startIdx = 1;
    for j =1:length(subjDirInfo)
       if ~isempty(strfind(subjDirInfo(j).name, 'anatomy'))
           subject(i).anatomy = 1; 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Restin_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Restin/rmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Restin(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Restin =unique(subject(i).Restin); 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Wrkmem_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Wrkmem/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Wrkmem(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Wrkmem =unique(subject(i).Wrkmem);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'StoryM_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/StoryM/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).StoryM(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).StoryM =unique(subject(i).StoryM);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Motort_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Motort/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Motort(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Motort =unique(subject(i).Motort);
       end
    end
end
clear tempDirInfo subjDirInfo subDirInfo sessionNum i j k;

workingDir = '/media/kuro/DATA1/MEGResult';