clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
roi_ids=cell2mat(roi_table.id);
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags_tested={-10:-4, 4:10,0 ,-10:10, -40:40};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase.mat' ],'b','r2','rnames');
        r2_perm=r2;
        b_perm=b;
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'b','r2','rnames');
        
        p=mean(r2_perm>r2,2);
         p(isnan(r2(:,1,1)))=NaN;
        pfdr=nan(size(p));
        [~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p)));
        pfdr=reshape(pfdr,size(p));
        
        pfwe=p*sum(~isnan(p(:)));
        
      
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','p','pfdr','r2_perm','b_perm','b');
        
        disp([exp ', lag' num2str(min(lags)) '-' num2str(max(lags))])
        rnames(pfdr<.05)
        nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05),  r2(pfdr<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_pfdr.nii']);
        
    end
end


