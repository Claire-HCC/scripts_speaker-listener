clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags_tested={-10:10,  -30:30};
permN=1000;
speakerSeed='vPCUN';


for ei=[3:4];%1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        
        load(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/SLg/perm/%s_lag%d-%d_permPhase',expdir,exp,timeUnit,speakerSeed,min(lags),max(lags)),'r','keptvox');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/' speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptT','keptvox');
        
        p=sum(r_perm>r,3)/permN;
        % p(regression_mask==0,:)=NaN;
        
        [peak peak_lagi]=max(r,[],2);
        peak_lags(:,1)=lags(peak_lagi);
        % fpcus on postive correlation
        p(peak<0)=NaN;
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        onsets_pfdr=nan(length(keptvox),1);
        onsets_p=nan(length(keptvox),1);
        for vi=1:length(keptvox);
            if sum(p_fdr(vi,:)<.05)>0;
                [lagi]=(find(p_fdr(vi,:)<.05));
                onsets_pfdr(vi)=lags(min(lagi));
            end
            
            if sum(p(vi,:)<.05)>0;
                [lagi]=(find(p(vi,:)<.05));
                onsets_p(vi)=lags(min(lagi));
            end
        end
        %        table(onsets(~isnan(onsets)),keptT(~isnan(onsets)))
        
        [~,~,onsetsOrder]=unique(onsets_pfdr);
        onsetsOrder(isnan(onsets_pfdr))=NaN;
        
        
        mat=nan(voxn,1);
        mat(keptvox)=peak_lags;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/'  speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak_lags.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=peak;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/'  speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=onsets_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/'  speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets_pfdr.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=onsets_p;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/' speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets_p.nii']);
        
                mat=nan(voxn,1);
                mat(keptvox)=onsetsOrder;
                nii=mat2nii(mat);
                save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/' speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsetsOrder_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/'  speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'keptT','lags','keptT','p','p_fdr','onsets_p','onsets_pfdr','peak','peak_lags','-v7.3');
        clear peak peak_lagi peak_lags
    end
end

