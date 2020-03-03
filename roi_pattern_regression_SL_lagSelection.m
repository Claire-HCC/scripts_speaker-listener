function roi_tr_pattern_regression_SL_lagSelection(perm)
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;

for ei=[1 2 4];%
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    listenerN=listenerNs(ei);
    rng(perm);
    [~,subj_perm]=sort(rand(listenerN,1));
    subj_test=subj_perm(1:(listenerN/3));
    subj_train=subj_perm(~ismember(subj_perm,subj_test));
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        for lagi=1:40;%
            lags=(-lagi):-1;
            keptT_s=find(([1:roi_tn]+min(lags))==1)+crop_start;
            keptT_e=min(roi_tn,roi_tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            g=nanmean(gdata(:,:,subj_train),3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            y=y-mean(y);
            
            g=nanmean(gdata(:,:,subj_test),3);
            y_test=g;
            y_test=y_test(:,keptT);
            y_test=y_test(:);
            y_test=y_test-mean(y_test);
            
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            X=X-mean(X);
            
            [b,~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
            
            r2_train(ri,li)=stats(1);
            r_test=y_test-[ones(size(X,1),1) X]*b;
            sst=sum((y_test-mean(y_test)).^2);
            ssr=sum(r_test.^2);
            r2_test(ri,li)=1-(ssr/sst);
            
            clear X
        end
    end
    
   save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/regression_SL_lagSelection_perm' num2str(perm)  ],'lags','rnames','keptT','r2_test','r2_train');
   clear r2_train r2_test
end

