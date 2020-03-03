% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
rname='vPCUN';
% -4 for merlin, -3 for sherlock
%
lags=-10:-1;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3%:4;%1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_s2l','r2_l','r2_ls','r2_s','eventLabels_LG_kept','rnames');
    eventLabels_LG_kept_sl=eventLabels_LG_kept;
    r2_s2l=F_s2l;%r2_ls-r2_l;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_event_lag0-0' ],'F','r2','eventLabels_LG_kept');
    eventLabels_LG_kept_ll=eventLabels_LG_kept;
    r2_ll=F;%nanmean(r2,4);
    
    % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    ris=find(ismember(rnames,rnames));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        
        % try dividing the story into early and late parts
        eventLabels_LG_kept=intersect(eventLabels_LG_kept_sl{ri},eventLabels_LG_kept_ll{ri});
        
        evs=(unique(eventLabels_LG_kept));
        % cols=distinguishable_colors(length(evs));
        cols=hot(length(evs)+10);
        
        r2_s2l_temp=r2_s2l{ri}(evs);
        r2_ll_temp=nanmean(r2_ll{ri}(evs),3);
        
        [r(ri,1) p(ri,1)]=corr(r2_s2l_temp,r2_ll_temp,'tail','right');
        
    end
    sigs_fdr=fdr0(p,0.05);
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) '_herd.mat'],'sigs_fdr','r','p');
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(sigs_fdr==1),'UniformOutput',0);
    roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
    nii=roiTable2wholeBrainNii_mor([roi_ids, r(sigs_fdr==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) '_fdr05.nii']);
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
