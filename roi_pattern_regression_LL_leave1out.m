function roi_pattern_regression_LL_leave1out(ei)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-15:15,-10:10,  -40:40, -10:-1 , -40:-1, 1:10,  1:40 };

% for ei=1;%2:4;
exp=experiments{ei};

for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata')
    [~ ,tn ,subjn]=size(gdata);
    
    keptT_s=min(find(([1:tn]+min(lags))>0))+crop_start;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([length(rnames) length(lags)+1 subjn]);
    r2=nan([length(rnames)  subjn]);
    F=nan([length(rnames)  subjn]);
    p=nan([length(rnames) subjn]);
    r2_byTime=nan([length(rnames),tn, subjn]);
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        disp(ri)
        if    exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            
            roi_voxn=size(gdata,1);
            
            for s=1:size(gdata,3)
                othersi=1:size(gdata,3);
                othersi=othersi(othersi~=s);
                others=nanmean(gdata(:,:,othersi),3);
                
                self=gdata(:,:,s);
                
                y=others(:,keptT);
                y=y(:);
                
                for li=1:length(lags);
                    X(:,:,li)=self(:,keptT+lags(li));
                end
                
                X=reshape(X,roi_voxn*length(keptT),length(lags));
                
                % centralized X
                X=X-mean(X);
                
                [b(ri,:,s),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                
                r2(ri,s)=stats(1);
                F(ri,s)=stats(2);
                p(ri,s)=stats(3);
                
                ymat=reshape(y,roi_voxn,length(keptT));
                rmat=reshape(r,roi_voxn,length(keptT));
                
                ssr=sum(rmat.^2);
                sst=sum((ymat-mean(y)).^2);
                r2_byTime(ri,:,s)=nan(tn,1);
                r2_byTime(ri,keptT,s)=1-ssr./sst;
                
                clear X
                
            end
        end
    end
    save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
    
    clear b F p r2 r2_byTime
end
% end


