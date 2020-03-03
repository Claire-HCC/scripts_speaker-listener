clear all;

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-10:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'rnames','keptT','r2_byTime_real','r2_byTime_null','r2_byTime_sig_fdr');
    r2_byTime_s=r2_byTime_real-mean(r2_byTime_null,3);
    sig_SL=r2_byTime_sig_fdr;
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag0-0_stats'],'r2_byTime_p','r2_byTime_sig_fdr','r2_byTime');
    listenerN=size(r2_byTime,3);
    r2_byTime_l=r2_byTime;
    r2_byTime_lm=mean(r2_byTime,3);
    sig_LL=r2_byTime_lm>0;
    
    for ri=9;
        figure
        hold on
        ti=(sig_SL(ri,:)==1 & sig_LL(ri,:)==1);
        scatter(r2_byTime_s(ri,ti)',r2_byTime_lm(ri,ti)',40,[1 0.7 0],'filled');
        
        ti=(sig_SL(ri,:)==1 & sig_LL(ri,:)==0);
        scatter(r2_byTime_s(ri,ti)',r2_byTime_lm(ri,ti)',40,[0 1 0],'filled');
        
        ti=(sig_SL(ri,:)==0 & sig_LL(ri,:)==0);
        scatter(r2_byTime_s(ri,ti)',r2_byTime_lm(ri,ti)',40,[0 0 1],'filled');
        
        ti=(sig_SL(ri,:)==0 & sig_LL(ri,:)==1);
        scatter(r2_byTime_s(ri,ti)',r2_byTime_lm(ri,ti)',40,[0.7 0.1 1],'filled');
        hold off
        
    end
     
end
