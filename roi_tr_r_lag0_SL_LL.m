
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
[~,rnames,~] = cellfun(@fileparts,cellstr(ls([expdir '/roi_mask/'  froidir '/mat/*.mat'])),'UniformOutput',0);

froidir='mor';
load([expdir '\roi_mask\' froidir '\roi_id_region.mat'],'roi_table');
tic % 15 min

load([expdir '/roi_mask/mor/roi_id_region.mat']);
for ei=2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        fl=[expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/simz_LL_' rname ];
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/simz_LL_' rname ],'sim','lags');
            simz_ll=sim;
            
            load([expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/simz_SL_' rname ],'sim','lags');
            simz_sl=sim;
            
            [r(ri,1) p(ri,1)]=corr(nanmean(simz_sl(lags==0,:,:),3)',nanmean( simz_ll(lags==0,:,:),3)');
            
            roi_id=cell2mat(roi_table.id(strmatch(rname,roi_table.region,'exact')));
            if ~isempty(roi_id)
                roi_ids(ri,1)=roi_id;
            else
                roi_ids(ri,1)=0;
            end
        end
    end
end

nii=roiTable2wholeBrainNii_mor([roi_ids r]);
save_nii(nii,[expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/r_lag0_SL_LL.nii']);

sig_pfwe=p<(0.05/sum(roi_ids>0));
nii=roiTable2wholeBrainNii_mor([roi_ids(sig_pfwe) r(sig_pfwe)]);
save_nii(nii,[expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/r_lag0_SL_LL_pfwe.nii']);

sig_pfdr=zeros(size(p));
sig_pfdr(p>0)=fdr0(p(p>0),0.05);
nii=roiTable2wholeBrainNii_mor([roi_ids(sig_pfdr) r(sig_pfdr)]);
save_nii(nii,[expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/r_lag0_SL_LL_pfdr.nii']);
