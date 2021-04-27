function roi_bined_pattern_regression_SL_leave1out(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[1 5 10 20 30]; % tr;
lags_tested={-10:10, -10:-3, -10:-1,-20:-3};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

% cropt start because there is clearly a speech-start effect in the
% listeners' data

% for ei=3:4;%1:2;
    exp=experiments{ei};
    mkdir(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/SL_leave1out_bined/perm/',expdir, exp, timeUnit, froidir));
    
    
    for binSizei=1:2;%length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=2;%1:length(lags_tested);
            lags=lags_tested{lagi};
            
            if ~exist(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/SL_leave1out_bined/binSize%d_lag%d-%d.mat',expdir, exp, timeUnit, froidir, binSize, min(lags),max(lags)));
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
                [~,tn,listenerN]=size(gdata);
                
                b=nan([length(rnames) tn length(lags)+1 listenerN]);
                r2=nan([length(rnames)  tn listenerN]);
                F=nan([length(rnames)  tn listenerN]);
                p=nan([length(rnames) tn listenerN]);
                
                for ri=1:length(rnames);
                    rname=rnames{ri};
                    
                    if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
                        gdata(:,:,subjects_excluded{ei})=NaN;
                        roi_voxn=size(gdata,1);
                        
                        for s=1:listenerN;
                            othersi=1:listenerN;
                            othersi=othersi(othersi~=s);
                            others=nanmean(gdata(:,:,othersi),3);
                            
                            for t=1:tn;
                                t_bin=t:(t+binSize-1);
                                
                                if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                                    % substract the global mean pattern
                                    y=others(:,t_bin);
                                    y=y(:);
                                    
                                    for li=1:length(lags);
                                        X(:,:,li)=data(:,t_bin+lags(li));
                                    end
                                    
                                    X=reshape(X,roi_voxn*length(t_bin),length(lags));
                                    
                                    % centralized X
                                    X=X-mean(X);
                                    
                                    % add an constant
                                    [b(ri,t,:,s),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                                    
                                    
                                    % https://www.riinu.me/2014/08/why-does-linear-model-without-an-intercept-forced-through-the-origin-have-a-higher-r-squared-value-calculated-by-r/
                                    r2(ri,t,s)=stats(1);
                                    F(ri,t,s)=stats(2);
                                    p(ri,t,s)=stats(3);
                                    
                                    clear X
                                end
                            end
                        end
                    end
                end
                save(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/SL_leave1out_bined/binSize%d_lag%d-%d',expdir, exp, timeUnit, froidir, binSize, min(lags),max(lags)),'b','F','r2','p','lags','rnames','binSize');
                
                clear b F p r2
            end
        end
    end
% end

beep
