function voxs2voxs_temporal_isc_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=20;
crop_end=20;

for ei=11;%1:12;%1:12;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/voxs2voxs//LL_leave1out/perm/']);
    
    load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/vox2vox//LL_leave1out/lag0-0'  ],'r','lags','keptvox','pfdr');
    rzm=squeeze(nanmean(atanh(r),3));
    % mask=(rzm>quantile(rzm(:),0.85));
    mask=(pfdr<.01) ;
    clear r
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata(:,:,subjects_excluded{ei})=NaN;
    gdata=gdata(mask==1,:,:);
    keptvox=keptvox(mask==1);
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    rzm=nan([length(keptvox)  length(keptvox)  ]);
    for si=1:listenerN;
        othersi=1:listenerN;
        othersi=othersi(othersi~=si);
        
        y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3);
        x=zscore(gdata(:,keptT,si),0,2);
        
        tmp = cat(3,rzm,atanh(corr(x',y')));
        rzm=nansum(tmp,3);
        
    end
    rzm=rzm./(listenerN-length(subjects_excluded{ei}));
    save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/voxs2voxs//LL_leave1out/isfc'  ],'rzm','lags','keptvox','keptT','-v7.3');
    
    [idx,C,sumd,D] = kmeans(rzm,6);
    
    mat=nan(voxn,1);
    mat(keptvox)=idx;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/voxs2voxs//LL_leave1out/networks_iscPfdrMask.nii']);
end



