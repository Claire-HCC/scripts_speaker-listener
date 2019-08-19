
clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

tic % 15 min
win_width=25;
win_step=1;
lags_new=-5:5;

for ei=1%:2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname ],'isc','lags');
            isc_sl=isc(ismember(lags,lags_new),:,:);
            load([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_LL_' rname ],'isc','lags');
            isc_ll=isc(ismember(lags,lags_new),:,:);
            
            load([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_SL_' rname  '_perm'],'isc_perm','lags');
            isc_sl_perm=isc_perm(ismember(lags,lags_new),:,:,:);
            load([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/iscz_LL_' rname '_perm'],'isc_perm','lags');
            isc_ll_perm=isc_perm(ismember(lags,lags_new),:,:,:);
            
            lags=lags_new;
            
            [isc_sl_bestT_r,isc_sl_bestT]=max(nanmean(nanmean(isc_sl,3),2));
            [isc_ll_bestT_r,isc_ll_bestT]=max(nanmean(nanmean(isc_ll,3),2));
            
            [isc_sl_perm_bestT_r,isc_sl_perm_bestT]=max(nanmean(nanmean(isc_sl_perm(:,:,:,:),3),2));
            [isc_ll_perm_bestT_r,isc_ll_perm_bestT]=max(nanmean(nanmean(isc_ll_perm(:,:,:,:),3),2));
            
            isc_sl_m=nanmean(isc_sl,3);
            isc_sl_perm_m=nanmean(isc_sl_perm,3);
            isc_ll_m=nanmean(isc_ll,3);
            isc_ll_perm_m=nanmean(isc_ll_perm,3);
            
            herd(ri,1)=corr(isc_sl_m(isc_sl_bestT,:)',isc_ll_m(lags==0,:)');
            peakT_sl(ri,1)=lags(isc_sl_bestT);
            
            iters = size(isc_sl_perm,4);
            for i=1:iters;
                isc_sl_perm_m_bestT(:,i)=isc_sl_perm_m(isc_sl_perm_bestT(i),:,:,i);
                isc_ll_perm_m_bestT(:,i)=isc_ll_perm_m(isc_ll_perm_bestT(i),:,:,i);
                
            end
            
            herd_perm(ri,:)=corr_col(isc_sl_perm_m_bestT,repmat(isc_ll_m(lags==0,:)',1,iters));
     
            
            roi_id=cell2mat(roi_table.id(strmatch(rname,roi_table.region,'exact')));
            if ~isempty(roi_id)
                roi_ids(ri,1)=roi_id;
            else
                roi_ids(ri,1)=0;
            end
            
        end
    end
     save([expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/herd_stats.mat'],'herd','herd_perm','roi_ids');
    %     nii=roiTable2wholeBrainNii_mor([roi_ids herd]);
    %     save_nii(nii,[expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/herd.nii']);
    %
    %        nii=roiTable2wholeBrainNii_mor([roi_ids peakT_sl]);
    %     save_nii(nii,[expdir '/' exp '/fmri/iscWin/' timeUnit '/roi/' froidir '/peakT_sl.nii']);
    
end


toc
