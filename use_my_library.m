function use_my_library(varargin)
if numel(varargin)<1
    error('Invalid option');
else
    library_name = lower(varargin{1});    
end
if numel(varargin)>1
    option = varargin{2};
    if ~(option == 1 ||option==0)
        error('Invalid option');
    end        
end
switch library_name
    case 'init'
        use_my_library('spm',0);
        use_my_library('ft',0);
    case 'spm'
        if option
            addpath('E:\matlabworks\spm12');
            addpath('E:\matlabworks\spm12\external\fieldtrip');
            spm('Defaults', 'eeg');
            warning('Use SPM -- Defualt is SPM12 --');
        else
            rmpath('E:\matlabworks\spm12');
            rmpath('E:\matlabworks\spm12\toolbox\DEM');
            clear ft_defaults
            rmpath('E:\matlabworks\spm12\external\fieldtrip',...
                'E:\matlabworks\spm12\external\fieldtrip\external\images',...
                'E:\matlabworks\spm12\external\fieldtrip\external\signal',...
                'E:\matlabworks\spm12\external\fieldtrip\connectivity',...
                'E:\matlabworks\spm12\external\fieldtrip\specest',...
                'E:\matlabworks\spm12\external\fieldtrip\plotting',...
                'E:\matlabworks\spm12\external\fieldtrip\inverse',...
                'E:\matlabworks\spm12\external\fieldtrip\forward',...
                'E:\matlabworks\spm12\external\fieldtrip\preproc',...
                'E:\matlabworks\spm12\external\fieldtrip\fileio',...
                'E:\matlabworks\spm12\external\fieldtrip\trialfun',...
                'E:\matlabworks\spm12\external\fieldtrip\statfun',...
                'E:\matlabworks\spm12\external\fieldtrip\external\fileexchange',...
                'E:\matlabworks\spm12\external\fieldtrip\utilities');
            rmpath('E:\matlabworks\spm12\external\bemcp',...
                'E:\matlabworks\spm12\external\ctf',...
                'E:\matlabworks\spm12\external\eeprobe',...
                'E:\matlabworks\spm12\external\mne',...
                'E:\matlabworks\spm12\external\yokogawa_meg_reader',...
                'E:\matlabworks\spm12\toolbox\dcm_meeg',...
                'E:\matlabworks\spm12\toolbox\spectral',...
                'E:\matlabworks\spm12\toolbox\Neural_Models',...
                'E:\matlabworks\spm12\toolbox\MEEGtools');
            warning('Remove SPM Path');
        end
    case 'ft'
        if option
            addpath('E:\matlabworks\fieldtrip\fieldtrip-20180918');
            ft_defaults
            warning('Use Fieldtrip -- Defualt is fieldtrip-20180918 --');
        else
            rmpath(genpath('E:\matlabworks\fieldtrip\fieldtrip-20180918'));
            warning('Remove Fieldtrip Path');
        end
    case 'fieldtrip'
        if option
            addpath('E:\matlabworks\fieldtrip\fieldtrip-20180918');
            ft_defaults
            warning('Use Fieldtrip -- Defualt is fieldtrip-20180918 --');
        else
            rmpath(genpath('E:\matlabworks\fieldtrip\fieldtrip-20180918'));
            warning('Remove Fieldtrip Path');
        end
    case 'bs'
        if option
            addpath(genpath('~\matlabwork\BrainStorm\brainstorm3'));
            warning('Use BrainStorm -- Defualt is BrainStorm3 --');
        else
            rmpath(genpath('~\matlabwork\BrainStorm\brainstorm3'));
            warning('Remove BrainStorm Path');
        end
    case 'mnet'
        if option
            addpath(genpath('E:\matlabworks\mnet\mnet0.92'));
            warning('Use MNET -- Defualt is mnet0.92 --');
        else
            rmpath(genpath('E:\matlabworks\mnet'));
            warning('Remove MNET Path');
        end
    case 'dodti'
        if option
            addpath(genpath('~\matlabwork\cdodti2.02_old'));
            warning('Use DoDTI -- Defualt is cDoDTI-2.02_old --');
        else
            rmpath(genpath('~\matlabwork\cdodti2.02_old'));
            warning('Remove DoDTI Path');
        end
    case 'dicm2nii'
        if option
            addpath(genpath('\projects1\dicm2nii'));
            warning('Use dicm2nii');
        else
            rmpath(genpath('\projects1\dicm2nii'));
            warning('Remove dicm2nii Path');
        end
    case 'ohba'
        if option
            addpath('\home\kuro\matlabwork\OHBA\osl-core');
            addpath('\home\kuro\matlabwork\OHBA\ohba-external');
            addpath('\home\kuro\matlabwork\OHBA\HMM-MAR');
            warning('Use OHBA toolbox');
        else
            rmpath('\home\kuro\matlabwork\OHBA\osl-core');
            rmpath('\home\kuro\matlabwork\OHBA\ohba-external');
            rmpath('\home\kuro\matlabwork\OHBA\HMM-MAR');
            warning('Remove OHBA toolbox');
        end
    case 'eeglab'
        if option
            addpath(genpath('E:\matlabworks\eeglab\eeglab14_1_2b'));
            warning('Use EEGLAB');
        else
            rmpath(genpath('E:\matlabworks\eeglab\eeglab14_1_2b'));
            warning('Remove EEGLAB Path');
        end
    otherwise
        warning('Invalid Option!')
end