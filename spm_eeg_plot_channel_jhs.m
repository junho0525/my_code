function spm_eeg_plot_channel_jhs(D, channellabel, trials);
% plot channel activity

% 2018-1-23 Junho Son
if ~numel(channellabel)==1
   channellabel = channellabel(1);
end
try
    data = D(:,:,:);
catch
end
ind = find(strcmp(D.chanlabels,channellabel));
data = squeeze(data(ind, :, trials));
time = D.time;
figure;
for t = trials
    plot(time, data(:,t)');
    hold on
end
hold off;