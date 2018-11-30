function D = spm_eeg_display_sensor_jhs(D)
% Diplay sensors SPM meeg Object in 2D coordinate

% 2018-1-23 Junho Son
try
    xylabel = D.chanlabels;
catch
    error('There is no ''chanlabels'' field in D');
end
try
    xy = coor2D(D);
catch
    error('Cannot read 2D coordinate of channels');
end

figure;
for n = 1:size(xy,2)
    text(xy(1,n), xy(2,n),xylabel(n));
end