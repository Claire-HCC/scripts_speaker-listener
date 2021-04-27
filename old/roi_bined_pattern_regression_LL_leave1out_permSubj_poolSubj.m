clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[10 20 30 40]; % tr;
lags_tested={0};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

for ei=1:2;
    exp=experiments{ei};
    
    for binSizei=1%:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1:length(lags_tested);
            lags=lags_tested{lagi};
            
            load(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/LL_leave1out_bined/binSize%d_lag0-0',expdir,exp,timeUnit,froidir,binSize),'r2');
            [roin,tn,listenerN]=size(r2);
            r2_llm=nan(size(r2));
            clear r2
            for perm=1:listenerN;
                if exist(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/LL_leave1out_bined/perm/binSize%d_lag0-0_permSL%03d.mat',expdir,exp,timeUnit,froidir,binSize,perm));
                    
                    load(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/LL_leave1out_bined/perm/binSize%d_lag0-0_permSL%03d',expdir,exp,timeUnit,froidir,binSize,perm),'r2','rnames');
                    r2_ll=r2;
                    r2_llm(:,:,perm)=nanmean(r2_ll,3);
                    
                end
            end
            r2=r2_llm;
            save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
            
        end
    end
end

