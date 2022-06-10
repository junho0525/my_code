function plot_meg_timeseries_jhs(data, options)
if isempty(data), error('No data!'); end
if nargin<2 || isempty(options), options.viewmode = 'butterfly'; options.ntrial = 5; options.trialoffset = 0; end
if ischar(data)
    try
        meg = load(data);
        if isfield(meg, 'data')
            datatype = 'ft';
            meg = meg.data;
        elseif isfield(meg, 'ftdata')
            datatype = 'ft';
            meg = meg.ftdata;
        elseif isfield(meg, 'D')
            datatype = 'spm';
            meg = spm_eeg_load(data);
        end
    catch
        use_my_library('ft',0);
        use_my_library('spm',1);
        datatype = 'spm';
        meg = spm_eeg_load(data);
    end
elseif isnumeric(data(:))
    datatype ='spm';
    meg = data;
elseif isfield(data, 'trial')
    datatype = 'ft';
    meg = data;
end

switch datatype 
    case 'ft'
        ntrial = min([length(meg.trial), options.ntrial, 5]);
        figure;
        for i = 1:ntrial
            subplot(ntrial, 1,i);
            plot(meg.time{i+options.trialoffset}, meg.trial{i+optioins.trialoffset});
        end
    case 'spm'
        ntrial = min([meg.ntrials, options.ntrial, 5]);
        figure;
        for i = 1:ntrial
            subplot(ntrial, 1,i);
            plot(meg.time, meg(:,:,i+options.trialoffset)');
        end
end

end