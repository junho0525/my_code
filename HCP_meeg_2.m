%mnet_hcp_meeg_defaults;
opengl hardware;
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'


% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

%experimentid = 'MEG';
%scanid = 'Resing';
subjectid = '100307';
experimentid=[subjectid '_MEG'];
mainpath=['/media/kuro/DATA1/HCPMEG/' experimentid];
megpath=fullfile(mainpath,[experimentid '_Restin_preproc'],subjectid,'MEG','Restin');
anatpath=fullfile(mainpath,[experimentid '_anatomy'],subjectid,'MEG','anatomy');
filename=fullfile(megpath,'rmegpreproc',[subjectid '_MEG_3-Restin_rmegpreproc.mat']);
resultprefix=fullfile(megpath,'icaclass',[subjectid '_MEG_3-Restin']);


cd(megpath)
hcp_read_matlab([resultprefix '_icaclass_vs.mat'])

% read the source and volume conduction model from current dir with
% outputs of previous pipelines



cd(anatpath);
hcp_read_matlab([subjectid '_MEG_anatomy_sourcemodel_2d.mat']);
sourcemodel2d=ft_convert_units(sourcemodel2d, 'cm');
sourcemodel2d.inside = 1:size(sourcemodel2d.pos,1);
sourcemodel2d.outside = [];
sourcemodelsubj = sourcemodel2d;

hcp_read_matlab(sprintf('%s.mat', [subjectid '_MEG_anatomy_headmodel']));
headmodel = ft_convert_units(headmodel, 'cm');
cd(megpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmeg=load(filename);
grad=rmeg.data.grad;
% gradBalanced = grad;
% gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
% grad=gradBalanced;
grad = ft_convert_units(grad, 'cm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the component data in order for ft_sourceanalysis to be able to
% swallow it
mixing   = comp_class.topo;
channels = comp_class.topolabel;
% normalisation of the topographies
for i = 1:size(mixing, 2)
  val(i) = 0.01*max(abs(mixing(:, i)));
  mixing(:, i) = mixing(:, i)/val(i);
end

% create a 'timelock' structure
tlck = [];
tlck.label = channels;
tlck.cov = eye(numel(tlck.label)); 
tlck.time=1;
tlck.grad = grad;
tlck.dimord = 'chan_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the forward solution

cfg = [];
cfg.vol = headmodel;
cfg.grid = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;

clear val tlck subjectid sourcemodelsubj sourcemodel2d resultprefix options mixing megpath;
clear rmeg mainpath i headmodel grad ft_default filename experimentid comp_class channels anatpath;

%%
cd ('~/Codes/M_code/Simul/HCP');
D = load('100307_MEG_3-Restin_rmegpreproc.mat');
D = spm_eeg_ft2spm(D.data, 'spm_meg');
S.D = D;
%D = spm_eeg_average(S);
clear S;
cd('..');
fid = fopen('100307_MEG_anatomy_transform.txt', 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;
D.inv = [];
D.inv{1}.forward = [];
D.inv{1}.method='Imaging';
D.inv{1}.forward.voltype='Single Shell';
D.inv{1}.forward.vol = ft_convert_units(cfg.vol,'m');
D.inv{1}.forward.modality = 'MEG';
D.inv{1}.forward.siunits=1;
D.inv{1}.forward.mesh=[];
D.inv{1}.forward.mesh_correction=[];
D.inv{1}.forward.mesh.face = cfg.grid.tri;
D.inv{1}.forward.mesh.vert = cfg.grid.pos;
D.inv{1}.forward.sensors = [];
D.inv{1}.forward.sensors = ft_convert_units(cfg.grad,'m');
D.inv{1}.forward.toMNI = transform.bti2spm;
D.inv{1}.forward.fromMNI = transform.spm2bti;
D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face;
D.inv{1}.mesh.tess_mni.vert = D.inv{1}.forward.mesh.vert;
Datareg(1).modality = 'MEG';
D.inv{1}.datareg(1) = Datareg(1);
clear Datareg;
D = spm_eeg_invert_ui(D);
%spm_eeg_invert_display(D, [-45,50,8]);