function roi_tr_bined_pattern_regression_LL

loc='cluster';
tic
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

% cropt start because ther eis clearly a spech-start effect in the
% listeners' data
binSize_tested=[10 15 30]; % tr;
lags=0;

for binSizei=1:2;%length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for ei=1:4;
        exp=experiments{ei};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
        tn=size(gdata,2);
        subjn=size(gdata,3);
        
        b=nan([length(rnames) tn length(lags)+1 subjn]);
        r2=nan([length(rnames)  tn subjn]);
        F=nan([length(rnames)  tn subjn]);
        p=nan([length(rnames) tn subjn]);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                roi_voxn=size(gdata,1);
                
                for s=1:subjn;
                    othersi=listeners(listeners~=s);
                    others=nanmean(gdata(:,:,othersi),3);
                    self=gdata(:,:,s);
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            y=others;
                            y=y(:,t_bin);
                            y=y(:);
                            
                            for li=1:length(lags);
                                X(:,:,li)=self(:,t_bin+lags(li));
                            end
                            
                            X=reshape(X,roi_voxn*length(t_bin),length(lags));
                            
                            % centralized X
                            X=X-mean(X);
                            
                            [b(ri,t,:,s),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                            
                            r2(ri,t,s)=stats(1);
                            F(ri,t,s)=stats(2);
                            p(ri,t,s)=stats(3);
                        end
                        clear X y
                    end
                end
                disp(ri)
            end
        end
        
        save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/LLselfother/regression_LL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','-v7.3');
        
        clear b F r2 p rnames
    end
end
toc





