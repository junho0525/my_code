function use_my_library(library_name, option)
library_name = lower(library_name);
if ~(option == 1 ||option==0)
    error('Invalid option');
end        
switch library_name
    case 'spm'
        if option
            addpath('/home/jhson/matlabwork/spm12');
            addpath('/home/jhson/matlabwork/spm12/external/fieldtrip');
            spm('Defaults', 'eeg');
            warning('Use SPM -- Defualt is SPM12 --');
        else
            rmpath('/home/jhson/matlabwork/spm12');
            rmpath('/home/jhson/matlabwork/spm12/toolbox/DEM');
            clear ft_defaults
            rmpath('/home/jhson/matlabwork/spm12/external/fieldtrip',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/external/images',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/external/signal',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/connectivity',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/specest',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/plotting',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/inverse',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/forward',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/preproc',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/fileio',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/trialfun',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/statfun',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/external/fileexchange',...
                '/home/jhson/matlabwork/spm12/external/fieldtrip/utilities');
            rmpath('/home/jhson/matlabwork/spm12/external/bemcp',...
                '/home/jhson/matlabwork/spm12/external/ctf',...
                '/home/jhson/matlabwork/spm12/external/eeprobe',...
                '/home/jhson/matlabwork/spm12/external/mne',...
                '/home/jhson/matlabwork/spm12/external/yokogawa_meg_reader',...
                '/home/jhson/matlabwork/spm12/toolbox/dcm_meeg',...
                '/home/jhson/matlabwork/spm12/toolbox/spectral',...
                '/home/jhson/matlabwork/spm12/toolbox/Neural_Models',...
                '/home/jhson/matlabwork/spm12/toolbox/MEEGtools');
            warning('Remove SPM Path');
        end
    case 'ft'
        if option
            addpath('/home/jhson/matlabwork/fieldtrip-20180206');
            addpath('/home/jhson/matlabwork/fieldtrip-20180206/template');
            warning('Use Fieldtrip -- Defualt is Fieldtrip-20180205 --');
        else
            rmpath('/home/jhson/matlabwork/fieldtrip-20180206');
            warning('Remove Fieldtrip Path');
        end
    case 'fieldtrip'
        if option
            addpath('/home/jhson/matlabwork/fieldtrip-20180206');
            addpath('/home/jhson/matlabwork/fieldtrip-20180206/template');
            warning('Use Fieldtrip -- Defualt is Fieldtrip-20180205 --');
        else
            rmpath('/home/jhson/matlabwork/fieldtrip-20180206');
            warning('Remove Fieldtrip Path');
        end
    case 'bs'
        if option
            addpath(genpath('~/matlabwork/BrainStorm/brainstorm3'));
            warning('Use BrainStorm -- Defualt is BrainStorm3 --');
        else
            rmpath(genpath('~/matlabwork/BrainStorm/brainstorm3'));
            warning('Remove BrainStorm Path');
        end
    case 'mnet'
        if option
            addpath('/home/jhson/matlabwork/mnet0.92');
            warning('Use MNET -- Defualt is mnet0.92 --');
        else
            rmpath(genpath('/home/jhson/matlabwork/mnet0.92'));
            warning('Remove MNET Path');
        end
    case 'dodti'
        if option
            addpath(genpath('/home/jhson/matlabwork/cdodti2.02'));
            warning('Use DoDTI -- Defualt is cDoDTI-2.02 --');
        else
            rmpath(genpath('/home/jhson/matlabwork/cdodti2.02'));
            warning('Remove DoDTI Path');
        end
    case 'dicm2nii'
        if option
            addpath(genpath('/projects1/dicm2nii'));
            warning('Use dicm2nii');
        else
            rmpath(genpath('/projects1/dicm2nii'));
            warning('Remove dicm2nii Path');
        end
    otherwise
        warning('Invalid Option!')
end