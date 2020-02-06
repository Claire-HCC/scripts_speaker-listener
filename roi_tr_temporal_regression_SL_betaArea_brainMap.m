clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=[1:4];;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification' ],'rnames','acc','sig_fdr_withinSigSL','b_s','b_l','p','p_adj','foldN');
        % sig_betaClass=sig_fdr_withinSigSL;
        
        blag_pos=nan([length(rnames) 1]);
        blag_neg=blag_pos;
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            % within roi showing signifcant difference in sl and ll beta
            % profiles.
            %    if sig_betaClass(ri)==1;
            b_s_temp=squeeze(b_s(ri,:,:))';
            b_l_temp=squeeze(b_l(ri,:,:))';
            
            % find the biggest peak within lag<=0
            [~,p,~,stats]=ttest(b_s_temp,0,'tail','right');
            
            b_norm=(nanmean(b_s_temp,1)/sum(nanmean(b_s_temp(:,lags<0 & p<.05),1)));
            a=lags.*b_norm.*(p<.05);
            blag_neg(ri,1)=nansum(a(lags<0));
            
             b_norm=(nanmean(b_s_temp,1)/sum(nanmean(b_s_temp(:,lags>0 & p<.05),1)));
            a=lags.*b_norm.*(p<.05);
            blag_pos(ri,1)=nansum(a(lags>0));
        end
    end
    blag_neg(blag_neg==0)=NaN;
    blag_pos(blag_pos==0)=NaN;
    save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea'  ],'rnames','blag_pos','blag_neg');
    
    nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( blag_neg)),   blag_neg(~isnan(  blag_neg))]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_neg.nii']);
    nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( blag_pos)),   blag_pos(~isnan(  blag_pos))]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_pos.nii']);
end


