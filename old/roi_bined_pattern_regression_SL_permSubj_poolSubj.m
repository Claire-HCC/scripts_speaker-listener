clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[1 5 10 20 30 40]; % tr;
lags_tested={-10:10 -10:-4, -10:-1, -20:-4, -30:-4,-10:-3};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

for ei=1:2;
    exp=experiments{ei};
    
    for binSizei=3%:2;%length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=[ 6];%1:length(lags_tested);
            lags=lags_tested{lagi};
            
            
            load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            r2_ll=r2;
            
            [roin,tn,listenerN]=size(r2_ll);
            
            r2_sl=nan(size(r2_ll));
            for perm=1:listenerN;
                if exist([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(perm) '.mat']);
                    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(perm)],'r2','rnames');
                    r2_sl(:,:,perm)=r2;
                    
                    %                     load(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/LL_leave1out_bined/perm/binSize%d_lag0-0_permSL%03d',expdir,exp,timeUnit,froidir,binSize,perm),'r2','rnames');
                    %                     r2_ll=r2;
                    %                     r2_llm=nanmean(r2_ll,3);
                    %                     herd_null(:,perm)=corr_col(r2_sl(:,keptT)',r2_llm(:,keptT)');
                end
            end
            r2=r2_sl;
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
            
        end
    end
end

