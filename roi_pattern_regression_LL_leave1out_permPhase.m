function roi_pattern_regression_LL_leave1out_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-15:15 -10:10,-4:4};
rname=rnames{ri};
permN=1000;

warning('off','stats:regress:RankDefDesignMat');

for ei=1:4;
    exp=experiments{ei};
    disp(ei)
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata')
        [~ ,tn ,subjn]=size(gdata);
        
        b=nan([1 length(lags)+1 subjn permN]);
        r2=nan([1  subjn permN]);
        F=nan([1  subjn permN]);
        p=nan([1 subjn permN]);
        r2_byTime=nan([1,tn, subjn, permN]);
        
        if  exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            
            [roi_voxn,tn,listenerN]=size(gdata);
            keptT_s=min(find(([1:tn]+min(lags))>0))+crop_start;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            for perm=1:permN;
                rng(perm);
                
                for si=1:listenerN;
                    if si==1;
                        [temp rp]=phase_rand2(gdata(:,:,si)',1);
                        gdata_perm(:,:,si)=temp';
                    else
                        gdata_perm(:,:,si)=(phase_rand2(gdata(:,:,si)',1,rp))';
                    end
                end
                
                for si=1:listenerN;
                    othersi=1:size(gdata,3);
                    othersi=othersi(othersi~=si);
                    others=nanmean(gdata_perm(:,:,othersi),3);
                    
                    self=gdata(:,:,si);
                    
                    y=others(:,keptT);
                    y=y(:);
                    
                    for li=1:length(lags);
                        X(:,:,li)=self(:,keptT+lags(li));
                    end
                    
                    X=reshape(X,roi_voxn*length(keptT),length(lags));
                    
                    % centralized X
                    X=X-mean(X);
                    
                    [b(1,:,si,perm),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                    
                    r2(1,si,perm)=stats(1);
                    F(1,si,perm)=stats(2);
                    p(1,si,perm)=stats(3);
                    
                    ymat=reshape(y,roi_voxn,length(keptT));
                    rmat=reshape(r,roi_voxn,length(keptT));
                    
                    ssr=sum(rmat.^2);
                    sst=sum((ymat-mean(y)).^2);
                    r2_byTime(1,:,si,perm)=nan(tn,1);
                    r2_byTime(1,keptT,si,perm)=1-ssr./sst;
                    
                    clear X
                    
                end
            end
        end
        save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname],'b','F','r2','p','lags','rnames','keptT','r2_byTime','-v7.3');
        
        clear b F p r2 r2_byTime
    end
end


