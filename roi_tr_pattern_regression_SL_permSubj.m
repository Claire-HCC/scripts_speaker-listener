function roi_tr_pattern_regression_SL_permSubj(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=1:4;
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
            
            roi_voxn=size(gdata,1);
            roi_tn=size(gdata,2);
            
            keptT_s=find(([1:roi_tn]+min(lags))==1)+crop_start;
            keptT_e=min(roi_tn,roi_tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            data_perm=gdata(:,:,perm);
            gdata_perm=gdata;
            gdata_perm(:,:,perm)=data;
            
            g=nanmean(gdata_perm,3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralize X
            X=X-mean(X);
            
            [b(ri,:),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
            
            r2(ri,1)=stats(1);
            F(ri,1)=stats(2);
            p(ri,1)=stats(3);
            
            ymat=reshape(y,roi_voxn,length(keptT));
            rmat=reshape(r,roi_voxn,length(keptT));
            ssr=sum(rmat.^2);
            sst=sum((ymat-mean(y)).^2);
            r2_byTime(ri,:)=nan(roi_tn,1);
            r2_byTime(ri,keptT)=1-ssr./sst;
            
            clear X
            disp(ri)
        end
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(perm)],'b','F','r2','r2_byTime','p','lags','rnames','keptT');
        clear b F p r2 r2_byTime
    end
end
toc
