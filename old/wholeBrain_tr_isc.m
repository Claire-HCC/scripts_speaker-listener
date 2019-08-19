
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;


tic % 15 min
win_width=25;
win_step=1;
lags=-20:20;

for ei=1%:2;
    exp=experiments{ei};
    
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker' ],'data','keptvox');
    
    
    speaker=data;
    
    
    for s=1:size(gdata,3);
        
        self=squeeze(gdata(:,:,s));
        othersi=1:size(gdata,3);
        othersi(othersi==s)=[];
        others=nanmean(gdata(:,:,othersi),3);
        
        isc_SL(:,1,s)=lagcorr_claire(self',speaker',lags);
        isc_LL(:,1,s)=corr_col(self',others');
    end
    
    isc=0.5*log((1+isc_SL)./(1-isc_SL));
    save([expdir '/' exp '/fmri/isc/' timeUnit '/wholeBrain/iscz_SL' ],'isc');
    mat=zeros(voxn,1);
    mat(keptvox)=nanmean(isc,3);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc/' timeUnit '/wholeBrain/iscz_SL.nii']);
    
    isc=0.5*log((1+isc_LL)./(1-isc_LL));
     save([expdir '/' exp '/fmri/isc/' timeUnit '/wholeBrain/iscz_LL' ],'isc');
     mat=zeros(voxn,1);
    mat(keptvox)=nanmean(isc,3);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc/' timeUnit '/wholeBrain/iscz_LL.nii']);
    
end
    
    toc
