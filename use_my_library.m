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
    case 'hcp'
        if option
            addpath(genpath('E:\matlabworks\megconnectome-3.0\megconnectome-3.0'));
            warning('Use megconnectome 3.0');
        else
            rmpath(genpath('E:\matlabworks\megconnectome-3.0\megconnectome-3.0'));
            warning('Remove megconnectome Path');
        end
    case 'spm'
        if option
            use_my_library('ft',0);
            addpath(genpath('E:\matlabworks\spm12'));
            addpath('E:\matlabworks\spm12\external\fieldtrip');
            spm('Defaults', 'eeg');
            warning('Use SPM -- Defualt is SPM12 --');
        else
            rmpath(genpath('E:\matlabworks\spm12'));
            warning('Remove SPM Path');
        end
     case 'spm-new'
        if option
            use_my_library('ft',0);
            addpath(genpath('E:\matlabworks\spm12-r7771'));
            addpath('E:\matlabworks\spm12-r7771\external\fieldtrip');
            spm('Defaults', 'eeg');
            warning('Use SPM12 --r7771 version--');
        else
            rmpath(genpath('E:\matlabworks\spm12-r7771'));
            warning('Remove SPM Path');
        end
    case 'ft'
        if option
            use_my_library('spm',0);
            addpath('E:\matlabworks\fieldtrip\fieldtrip-20201128');
            ft_defaults
            warning('Use Fieldtrip -- Defualt is fieldtrip-20201128 --');
        else
            rmpath(genpath('E:\matlabworks\fieldtrip\fieldtrip-20201128'));
            warning('Remove Fieldtrip Path');
        end
    case 'fieldtrip'
        if option
            use_my_library('spm',0);
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
            addpath('E:\matlabworks\mnet\mnet0.92');
            addpath('E:\matlabworks\mnet\mnet0.92\labels');
            addpath('E:\matlabworks\mnet\mnet0.92\atlas');
            addpath('E:\matlabworks\mnet\mnet0.92\codes')
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
    case 'hmmmar'
        if option
            addpath(genpath('E:\matlabworks\osl\HMM-MAR-master'));
            addpath(genpath('E:\matlabworks\osl\osl-core-master'));
            addpath(genpath('E:\matlabworks\osl\ohba-external-master'));
            warning('Use HMM-MAR toolbox');
        else
            rmpath(genpath('E:\matlabworks\osl\HMM-MAR-master'));
            rmpath(genpath('E:\matlabworks\osl\osl-core-master'));
            rmpath(genpath('E:\matlabworks\osl\ohba-external-master'));
            warning('Remove HMM-MAR toolbox');
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