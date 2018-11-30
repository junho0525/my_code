%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descriptions have subjectId and sessionId for DCM for each DCMs cell    %
% lines.                                                                  %
% A pair of DCMs that have same subjectId and sessionId inserted in       %
% a single cell line.                                                     %
% The first one is DCM_Alpha_High, and the second one is DCM_Alpha_Low.   %
% You should specify worikingDir for your own system.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Subject list from HCP Dataset
workingDir = '.';
dirInfo = dir(workingDir);
fileRemoveIdx = zeros(0,1);
for i =1:length(dirInfo)
    fileRemoveIdx(end+1) = isempty(strfind(dirInfo(i).name, 'DCM'));
end
dirInfo(find(fileRemoveIdx)) = []; % get filenames contain string 'DCM'.
clear fileRemoveIdx i

Descriptions = [];                                                         
DCMs = cell(0,2);                                                          

i = 1;
while i < length(dirInfo)
    fileNameSplit1 = strsplit(dirInfo(i).name,'_');
    fileNameSplit2 = strsplit(dirInfo(i+1).name,'_');
    if strcmp(fileNameSplit1{2},fileNameSplit2{2}) && strcmp(fileNameSplit1{3},fileNameSplit2{3})
        Descriptions(end+1).subjectId = fileNameSplit1{2};
        Descriptions(end).sessionId = fileNameSplit1{3};
        load(fullfile(workingDir,dirInfo(i).name));
        load(fullfile(workingDir,dirInfo(i+1).name));
        DCMs(end+1,:) = {DCM_High_Alpha,DCM_Low_Alpha};
    else
        i=i+1;
        continue;
    end
    i = i+2;
end

%% Sort DCMs By sessionId
Descriptions_new= [];                                                      
i = 1;
while i < length(Descriptions)
    if~(strcmp(Descriptions(i).subjectId, Descriptions(i+1).subjectId)&&...
            strcmp(Descriptions(i+1).subjectId,Descriptions(i+2).subjectId))
        i=i+1;
        continue;
    end
    Descriptions_new(end+1).subjectId =Descriptions(i).subjectId;
    i = i+3;
end

