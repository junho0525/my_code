function sourceout = mnet_source_reconstruction(cfg, data, sourcemodelfile, headmodelfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG Source Localization Code                                            %
%     cfg.method          - specify source localize method                %
%         time domain analysis
%             'lcmv', 'eloreta','sloreta'
%         frequency domain analysis
%             'pcc', 'dics'
%
%     cfg.resample_trials - specify method for resampling trials          %
%             'no'        - average the data if multiple trials are       %
%                           present and no resampling is neccessary       %
%         'permutation'   - compute the average and covariance with       %
%                           random permutation                            %
%        'randomization'  - compute the average and covariance with       %
%                           random resampling                             %
%          'jackknife'    - compute the jackknife repetitions for         %
%                           the average and covariance                    %
%          'bootstrap'    - compute the bootstrap repetitions for         %
%                           the average and covariance                    %
%          'rawtrial'     - do not do any resampling, keep the            %
%                           single-trial covariances and single-trial ERFs%
%                           (rename them to avg for convenience)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.09.06 16:22 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Leadfield
load(headmodelfile);
load(sourcemodelfile);

tmpcfg = [];
tmpcfg.headmodel = headmodel;
tmpcfg.grid = sourcemodel2d;
tmpcfg.grad = data.grad;
tmpcfg.channel = data.label;

leadfield = ft_prepare_leadfield(tmpcfg);
%% Source Estimation.
tmpcfg                   = [];
tmpcfg.grid              = leadfield;
tmpcfg.headmodel         = headmodel;
switch cfg.resample_trials
    case 'no'
        
    case 'permutation'
        tmpcfg.permutation = 'yes';
    case 'randomization'
        tmpcfg.randomization = 'yes';
    case 'jackknife'
        tmpcfg.jackknife = 'yes';
    case 'bootstrap'
        tmpcfg.bootstrap = 'yes';
    case 'singletrial'
        error('fieldtrip-20180208 does not support this option yes');
    case 'rawtrial'
        tmpcfg.rawtrial = 'yes';
    otherwise
        error('Undefined option on cfg.resample_trials');
end
switch cfg.method
    case 'lcmv'
        tmpcfg.method = 'lcmv';
        tmpcfg.lcmv.lambda = '10%';
        tmpcfg.lcmv.fixedori = 'yes';
        tmpcfg.lcmv.projectnoise = 'yes';
        tmpcfg.lcmv.keepleadfield = 'yes';
        source = ft_sourceanalysis(tmpcfg, data);
        sourceout = source;
        source = ft_sourcedescriptives([], source);
    case 'eloreta'
        error('This option is not implemented at mnet 0.92 - 20180906');
    case 'sloreta'
        error('This option is not implemented at mnet 0.92 - 20180906');
    case 'pcc'
        error('This option is not implemented at mnet 0.92 - 20180906');
    case 'dics'
        error('This option is not implemented at mnet 0.92 - 20180906');
    otherwise
        error('Undefined option on cfg.method');
end

%% Visualize Source Activity
figure;
bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;
subplot(2,2,1);
for trial = 1:length(data.trial)
    if ~any(any(isnan(data.trial{trial})))&&~isempty(data.trial{trial})
        plot(data.time{trial},data.trial{trial});
        title(['Channel Signal of Trial ' num2str(trial)]);
        break;
    end
end


switch cfg.method
    case 'lcmv'
        subplot(2,2,3);
        for i = 1:length(sourceout.trial(trial).mom)
            if ~any(isnan(sourceout.trial(trial).mom{i}))&&~isempty(sourceout.trial(trial).mom{i})
                plot(source.time,sourceout.trial(trial).mom{i});
                hold on
            end
        end
        hold off
        title(['Source Signal of Trial ' num2str(trial) ' (Beamformer LCMV)']);
        subplot(2,2,[2 4]);
        scale_factor = median(floor(log10(source.avg.nai(find(~isnan(source.avg.nai))))));
        source.avg.nai = 10.^(log10(source.avg.nai)-scale_factor);
        ft_plot_mesh(bnd);
        ft_plot_mesh(bnd, 'vertexcolor', source.avg.nai, 'facealpha',source.avg.nai,'alphalim','auto','colormap','jet');
        lighting gouraud
        camlight
        title(sprintf('Source Reconstruction of Data\n(Beamformer LCMV)'));
end


