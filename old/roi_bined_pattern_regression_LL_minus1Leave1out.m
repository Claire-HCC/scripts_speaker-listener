function roi_bined_pattern_regression_LL_minus1Leave1out(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

% cropt start because ther eis clearly a spech-start effect in the
% listeners' data
binSize_tested=[1 5 10 20 30 ]; % tr;
lags=0;

exp=experiments{ei};
for binSizei=1:2;%length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    if ~exist([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ]);
        
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    
    b=nan([length(rnames) tn length(lags)+1 listenerN listenerN]);
    r2=nan([length(rnames)  tn listenerN listenerN]);
    F=nan([length(rnames)  tn listenerN listenerN]);
    p=nan([length(rnames) tn listenerN listenerN]);
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            roi_voxn=size(gdata,1);
            gdata(:,:,subjects_excluded{ei})=NaN;
            
            for fakeSpeaker=1:listenerN;
                for s=1:listenerN;
                    if s~=fakeSpeaker;
                        othersi=1:listenerN;
                        othersi=othersi(~ismember(othersi,[s fakeSpeaker]));
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
                                
                                [b(ri,t,:,fakeSpeaker,s),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                                
                                % https://www.riinu.me/2014/08/why-does-linear-model-without-an-intercept-forced-through-the-origin-have-a-higher-r-squared-value-calculated-by-r/
                                r2(ri,t,fakeSpeaker,s)=stats(1);%1-sum(r.^2)/sum((y).^2);
                                F(ri,t,fakeSpeaker,s)=stats(2);
                                p(ri,t,fakeSpeaker,s)=stats(3);
                                
                                clear X
                            end
                        end
                    end
                end
            end
        end
    end
    mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/']);
    save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','-v7.3');
    
    clear b F r2 p
    end
end







