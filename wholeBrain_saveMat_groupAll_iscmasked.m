loc='mypc';
set_parameters;
timeUnit='tr' ;
perc=0.30;
load([expdir 'roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain_gm.mat'],'roimask')

for ei=5;%1:11;%
    exp=exp_parameters.experiments{ei};
    
    %  if exist( sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/isc%dPercMask.nii',expdir,exp,timeUnit,perc*100))==0;
    clear gdata
    
    % within the gray matter mask, keep voxels showing top 30% isc
    if ismember(exp,'pieman_rest');
        load([expdir '/' 'pieman_old' '/fmri/temporal/circularlagcorr/' timeUnit '/vox/LL_leave1out/isc'   ],'rzm','lags','keptvox');
    else
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox/LL_leave1out/isc'   ],'rzm','lags','keptvox');
    end
    isc=rzm(:,lags==0);
    % keep the top 30% voxels
    thr=quantile(isc,1-perc);
    mask=zeros(voxn,1);
    mask(intersect(keptvox(isc>thr) , find(roimask==1)))=1;
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata_listener=gdata;
    keptvox_wholeBrain_L=keptvox;
    
    gdata=gdata_listener(logical(mask(keptvox_wholeBrain_L)) ,:,:);
    keptvox=keptvox_wholeBrain_L(logical(mask(keptvox_wholeBrain_L)));
    
    fout= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,exp,timeUnit,perc*100);
    save(fout,'gdata','keptvox','mask','-v7.3');
    
    %      nii=mat2nii(mask);
    %     save_nii(nii, sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/isc%dPercMask.nii',expdir,exp,timeUnit,perc*100));
    % end
end
