function roi_tr_pattern_lagcorr_SLeach_betaClassification(ei)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};
foldN=2;
bootsSampleN=10000;
exp=experiments{ei};
listeners_exp=listeners{ei};

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/lagcorr_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats' ],'sig_fdr');
    sig_SL=sig_fdr;
    
    load(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/SLeach/perm/lagcorr_SLeach_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),0),'r');
    
    b_s=b(:,2:end,:);
    subjN=size(b_s,3);
    subsetSize=subjN/foldN;
    
    b_l=nan([size(b_s) subjN]);
    for si=1:length(listeners_exp);
        perm=listeners_exp(si);
        load(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/SLeach/perm/lagcorr_SLeach_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b');
        b_l(:,:,si,:)=reshape(b(:,2:end,:),[size(b,1) length(lags) 1 subjN]);
    end
    
    for booti=1:bootsSampleN;
        subj_shuffled=randperm(subjN);
        subj_test=subj_shuffled(1:(subjN/foldN));
        subj_train=find(~ismember(1:subjN,subj_test));
        
        b_s_test=b_s(:,:,subj_test);
        b_s_train=b_s(:,:,subj_train);
        
        b_l_test=b_l(:,:,subj_test,subj_test);
        b_l_test=reshape(b_l_test,size(b_l_test,1),size(b_l_test,2),length(subj_test)^2);
        b_l_train=b_l(:,:,subj_train,subj_train);
        b_l_train=reshape(b_l_train,size(b_l_train,1),size(b_s_train,2),length(subj_train)^2);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]) ;
                
                speaker_train=squeeze(nanmean(b_s_train(ri,:,:),3));
                listener_train=squeeze(nanmean(b_l_train(ri,:,:),3));
                
                data_test=[squeeze(b_s_test(ri,:,:))' ;squeeze(b_l_test(ri,:,:))' ];
                label_test=[zeros(size(b_s_test,3),1); ones(size(b_l_test,3),1)];
                
                label_predicted=corr_col(data_test',repmat(listener_train',1,size(data_test,1)))>corr_col(data_test',repmat(speaker_train',1,size(data_test,1)));
                acc(ri,booti)=mean(label_test==label_predicted');
            else
                acc(ri,booti)=NaN;
            end
        end
    end
    load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/lagcorr_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification' ],'rnames','acc','b_s','b_l','foldN');
    
    
      nii=roiTable2wholeBrainNii_mor([roi_ids  mean(acc,2)]);
    save_nii(nii,['Y:\claire\speaker-listener\' exp '\fmri\pattern_regression\tr\roi\mor\SLeach\regression_SLeach_lag-10-10_classification_acc.nii'])

    
    p=nansum(acc<0.5,2)/bootsSampleN;
    ris=find(sum(isnan(acc),2)==0);
    p_adj=nan(size(p));
    [~,~,p_adj(ris)]=fdr(p(ris));
    
    ris=find(sum(isnan(acc),2)==0 & sig_SL==1);
    p_adj_withinSigSL=nan(size(p));
    [~,~,p_adj_withinSigSL(ris)]=fdr(p(ris));
    sig_fdr_withinSigSL=(p_adj_withinSigSL<.05);
    rnames(sig_fdr_withinSigSL==1)
    
    save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/lagcorr_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification' ],'rnames','acc','sig_fdr_withinSigSL','b_s','b_l','p','p_adj','foldN','p_adj_withinSigSL');
end
