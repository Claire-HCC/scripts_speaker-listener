function roi_pattern_regression_SL_g(ei)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:-4. 4:10,0 ,-10:10, -40:40};


exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/perm']);

for lagi=4;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rnames{1} '.mat' ],'data')
    tn=size(data,2);
    
    b=nan([length(rnames) length(lags)+1 ]);
    r2=nan([length(rnames) 1 ]);
    F=nan([length(rnames) 1]);
    p=nan([length(rnames) 1]);
    r2_byTime=nan([length(rnames),tn]);
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            gdata(:,:,subjects_excluded{ei})=NaN;
            roi_voxn=size(gdata,1);
            
            keptT_s= max(-min(lags),0)+crop_start+1;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            g=nanmean(gdata,3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            X=X-mean(X);
            
            [b(ri,:),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
            
            r2(ri,1)=stats(1);
            F(ri,1)=stats(2);
            p(ri,1)=stats(3);
            
            ymat=reshape(y,roi_voxn,length(keptT));
            rmat=reshape(r,roi_voxn,length(keptT));
            ssr=sum(rmat.^2);
            sst=sum((ymat-mean(y)).^2);
            
            r2_byTime(ri,:)=nan(tn,1);
            r2_byTime(ri,keptT)=1-ssr./sst;
            
            clear X
        end
    end
    
    save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
    clear b F p r2 r2_byTime
end


