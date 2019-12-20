% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
% -4 for merlin, -3 for sherlock
lags_tested={-10:-1,-10:-4};
binSize_tested=[  10 15 30 ];
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

% load([expdir '/scripts_speaker-listener/merlin_' rname '_hmm_findListenersEventInSpeaker.mat'], 'segments_LG');

%r2_l leads to lower r_ls-r2_l; which might not be fair to the speaker,
%because listner can predict, but not as good as speaker.
% the solutions are:
% 1. set lags=-10:-4
% 2. use F-value as coupling index
% don't use granger causality.

for ei=1:4;
    exp=experiments{ei};
    fsize=[20 9];
    % f1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    % f2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    
    for bi=1:length(binSize_tested);
        binSize=binSize_tested(bi);
        binStep=binSize;
        
        for lagi=1:length(lags_tested);
            lags=lags_tested{lagi};
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_s2l','F_l','r2_l','r2_ls','rnames');
            coupling_sl=r2_ls-r2_l;
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'F','r2');
            coupling_ll=nanmean(r2,3);
            listenerN=size(r2,3);
            
            keptT=find(sum(isnan(coupling_sl),1)==0  & sum(isnan(coupling_ll),1)==0 );
            keptT=keptT(1:binStep:end);
            
            herd(:,1)=corr_col(coupling_sl(:,keptT)',coupling_ll(:,keptT)');
            
            for s=1:listenerN;
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(s)],'F_s2l','r2_l','r2_ls','rnames');
                coupling_sl_perm=r2_ls-r2_l;
                
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag0-0_perm' num2str(s)],'F','r2');
                coupling_ll_perm=nanmean(r2,3);
                
                herd_null(:,s)=corr_col(coupling_sl_perm(:,keptT)',coupling_ll_perm(:,keptT)');
                
            end
            herd_null_z=atanh(herd_null);
            herd_z=atanh(herd);
            for ri=1:length(rnames);
                [~,p(ri,1)]=ttest(herd_null_z(ri,:),herd_z(ri),'tail','left');
            end
            
            sig_fdr=fdr0(p,0.05);
           save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'p','herd_null','herd','rnames','lags','keptT','sig_fdr');
            
            rnames(sig_fdr==1)
            
            figure;
            boxplot(herd_null');
            % scatter(repmat([1:length(rnames)],1,listenerN),herd_null(:),40,'k','filled');
            hold on
            scatter(1:length(rnames),herd,40,'r','filled');
            hold off
            set(gca,'xticklabels',rnames,'fontsize',8)
            xtickangle(70)
            grid on
            title([exp ',binsize' num2str(binSize) ', lag' num2str(min(lags)) '-' num2str(max(lags))])
        end
    end
end
