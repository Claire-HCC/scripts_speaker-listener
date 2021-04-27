function roi_temporal_lagcorr_classification(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};
foldN=2;
bootsSampleN=10000;
exp=experiments{ei};

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    %  load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats' ],'sig_fdr');
    %  sig_SL=sig_fdr;
    
    load(sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/SL_each/lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'r');
    rz_s=atanh(r);
    rz_s(:,:,subjects_excluded{ei})=NaN;
    
    load(sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/LL_each/lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'r');
    rz_l=atanh(r);
    
    [~,~,listenerN]=size(rz_s);
    
    for booti=1:bootsSampleN;
        subj_shuffled=randperm(listenerN);
        subj_test=subj_shuffled(1:(listenerN/foldN));
        subj_train=find(~ismember(1:listenerN,subj_test));
        
        rz_s_test=rz_s(:,:,subj_test);
        rz_s_train=rz_s(:,:,subj_train);
        
        rz_l_test=rz_l(:,:,subj_test,subj_test);
        rz_l_test=reshape(rz_l_test,size(rz_l_test,1),size(rz_l_test,2),length(subj_test)^2);
        rz_l_train=rz_l(:,:,subj_train,subj_train);
        rz_l_train=reshape(rz_l_train,size(rz_l_train,1),size(rz_s_train,2),length(subj_train)^2);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]) ;
                
                speaker_train=squeeze(nanmean(rz_s_train(ri,:,:),3));
                listener_train=squeeze(nanmean(rz_l_train(ri,:,:),3));
                
                data_test=[squeeze(rz_s_test(ri,:,:))' ;squeeze(rz_l_test(ri,:,:))' ];
                label_test=[zeros(size(rz_s_test,3),1); ones(size(rz_l_test,3),1)];
                
                label_predicted=corr_col(data_test',repmat(listener_train',1,size(data_test,1)))>corr_col(data_test',repmat(speaker_train',1,size(data_test,1)));
                acc(ri,booti)=mean(label_test==label_predicted');
            else
                acc(ri,booti)=NaN;
            end
        end
    end
    
    p=nanmean(acc<.5,2);
    pfdr=nan(size(p));
    [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
    rnames(pfdr<.05)
    save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification' ],'rnames','acc','rz_s','rz_l','p','pfdr','foldN');
end
