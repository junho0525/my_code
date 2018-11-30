function [source, weights, leadfield,cfg] = bf_inverse_lcmv_jhs(data, config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG/EEG Source Localization Using Beamformer LCMV                       %
% Original code is from HCP_beamform.m code in OHBA toolbox               %
%                                                                         %
% Essential arguments                                                     %
%    Input Data    - fieldtrip data format                                %
%    Source Model  - fieldtrip source model                               %
%    Head Model    - fieldtrip Volume conductance model                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.07.03 11:23 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if strcmp(ft_datatype(data),'unknown'),error('Data is not a fieldtrip datatype!');end
if ~isfield(config, 'sourcemodel'),error('There is no sourcemodel!');end
if ~isfield(config, 'headmodel'),error('There is no headmodel!');end
%% 
% Compute the leadfields
cfg             = [];
cfg.vol         = config.headmodel;
cfg.grid        = config.sourcemodel;
cfg.grad        = data.grad;
cfg.channel     = data.label;
cfg.normalize   = 'no';
cfg.reducerank  = 2;

% prepare leadfields
[leadfield,cfg] = ft_prepare_leadfield(cfg);

% Compute the beamformer weights
C = 0;
for trl = 1:numel(data.trial)
    C  = C + osl_cov(data.trial{trl});    
end
C = C ./ numel(data.trial);

% Estimate rank of data covariance matrix:
eigDiff = diff(log(svd(C)));
rankVec = eigDiff<10*median(eigDiff);
rankVec(1:200) = false;
rankC = min([find(rankVec,1,'first'),rank(C)]);

invC = pinv_plus(C,rankC); 

% rank for leadfields
rankLF = 2;


% Compute weights for each voxel
weights = zeros(size(leadfield.leadfield,2), size(leadfield.leadfield{1},1));
emptyarray = [];
nanarray = [];
for voxel = 1:length(weights)
    lf  = leadfield.leadfield{voxel};
    if isempty(lf) 
        emptyarray(end+1) = voxel;
        continue;
    end
    if isnan(lf)
        nanarray(end+1) = voxel;
        continue;
    end
    % Scalar beamformer - reduce leadfield to rank 1
    [u, ~] = svd(real(pinv_plus(lf' * invC *lf, rankLF, 0)),'econ');
    eta = u(:,1);
    lf = lf * eta;
    
    % LCMV weights
    weights(voxel,:) = pinv_plus(lf' * invC * lf, rankLF, 0) * lf' * invC;
end
source = zeros(size(leadfield.leadfield,2),size(data.trial{1},2), length(data.trial));
for trial = 1:length(data.trial)
    source(:,:,trial) = weights*data.trial{trial};
end
end