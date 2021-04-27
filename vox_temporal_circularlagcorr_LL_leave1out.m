function vox_temporal_circularlagcorr_LL_leave1out(ei)

loc='mypc';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;

exp=exp_parameters.experiments{ei};

mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox//LL_leave1out/perm/']);

f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
load(f,'gdata','keptvox');
gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;

[~,tn,listenerN]=size(gdata);
keptT=(crop_start+1):(tn-crop_end);
lags=0;%-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);

r=nan([length(keptvox)  length(lags) listenerN   ]);

for si=1:listenerN;
    othersi=1:listenerN;
    othersi=othersi(othersi~=si);
    
    y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3)';
    x=zscore(gdata(:,keptT,si),0,2)';
    r(:,:,si)=circularlagcorr_byCol(x,y,lags)';
    disp(si)
end
rzm=nanmean(atanh(r),3);

save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox/LL_leave1out/isc'   ],'rzm','lags','keptvox','keptT','-v7.3');

mat=nan(voxn,1);
mat(keptvox,1)=atanh(rzm(:,lags==0));
nii=mat2nii(mat);
save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox/LL_leave1out/isc_r.nii'   ])