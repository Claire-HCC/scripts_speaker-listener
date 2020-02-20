
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;
speakerSeed='vPCUN';

for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        
        load(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/%s/SLg/perm/%s_lag%d-%d_permPhase',expdir,exp,timeUnit,froidir,speakerSeed,min(lags),max(lags)),'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/' froidir '/SLg/' speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptT','keptvox');
        p=sum(r_perm>r,3)/permN;
        %  p(regression_mask==0,:)=NaN;
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        p_fdr=p;
        onsets=nan(length(keptvox),1);
        onsets1=nan(length(keptvox),1);
        onsets2=nan(length(keptvox),1);
        for vi=1:length(keptvox);
            
            if sum(p_fdr(vi,:)<.05)>0;
                [lagi]=(find(p_fdr(vi,:)<.05));
                onsets(vi)=lags(min(lagi));
                
                try
                    onsets1(vi)=min(intersect(lags(lagi),lags(lags<0)));end
                try
                    onsets2(vi)=min(intersect(lags(lagi),lags(lags>=0)));end
            end
        end
        % table(onsets(~isnan(onsets)),keptT(~isnan(onsets)))
        
        %         [~,~,onsetsOrder]=unique(onsets);
        %         onsetsOrder(isnan(onsets))=NaN;
        
        mat=nan(voxn,1);
        mat(keptvox)=onsets;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets.nii']);
        %
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsetsOrder)),   onsetsOrder(~isnan( onsetsOrder))]);
        %         save_nii(nii,[expdir '/' exp '/fmvi/temporal_lagcorr/' timeUnit '/roi2wholeBrain/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsetsOrder.nii']);
        %
        %         save([expdir '/' exp '/fmvi/temporal_lagcorr/' timeUnit '/roi2wholeBrain/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'keptT','r','lags','keptT','r_perm','p','p_fdr','onsets','onsetsOrder');
    end
end

