clear all
loc='mypc';

set_parameters;
perc=0.15;
experiments={'sherlock'    'merlin'    'pieman_old'    'black'    'forgot'    'ABC'    'bronx'    'pieman'};

for ei=1:length(experiments);
    exp=experiments{ei};
    
    % masked by story isc
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,exp,'tr',perc*100);
    load(f,'mask');
    if ei==1;
        masks=mask;
    else
    masks=masks+mask;
    end
end
masks=masks/length(experiments);
nii=mat2nii(masks);
save_nii(nii,sprintf('%s/roi_mask/isc%dPercMasks.nii',expdir,perc*100));

roimask=masks;
save(sprintf('%s/roi_mask/isc%dPercMasks.mat',expdir,perc*100),'roimask','experiments');

