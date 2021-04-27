function vox2vox_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
lags_tested={-20:20 };

for ei=[ 1 2 4 9:12]%:12;%1:12;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/vox2vox//LL_leave1out/perm/']);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata(:,:,subjects_excluded{ei})=NaN;
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        r=nan([length(keptvox)  length(lags) listenerN ]);
        
        for si=1:listenerN;
            othersi=1:listenerN;
            othersi=othersi(othersi~=si);
            
            y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3);
            x=zscore(gdata(:,keptT,si),0,2);
            
            r(:,:,si)=(lagcorr_claire(x',y',lags))';
            
        end
        
        rz=atanh(squeeze(r));
        t=[];
        p=[];
        for lagi=1:length(lags);
            [~,p(:,lagi),~,stats]=ttest(squeeze(rz(:,lagi,:))');
            t(:,lagi)=stats.tstat;
        end
        [~,~,pfdr]=fdr(p(:));
        pfdr=reshape(pfdr,size(p));
        
        rzm=squeeze(nanmean(rz,3));
        [peaks, peakLagi]=max(rzm');
        
        peakLags=lags(peakLagi);
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/vox2vox//LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags))  ],'r','lags','keptvox','keptT','t','p','pfdr','rzm','peaks','peakLags','-v7.3');
        
        %              peakLag0_mask=zeros(voxn,1);
        %     peakLag0_mask(keptvox)=(peakLags==0 & peaks>0);
        %   nii=mat2nii( peakLag0_mask);
        %  save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/vox2vox//LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) '_peakLag0_mask.nii'])
        
    end
end


