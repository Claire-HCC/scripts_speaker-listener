function roi_pattern_regression_SL_each(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-40:40,-10:-4. 4:10, -10:10, -10:-1, 0, 1:10};
exp=experiments{ei};
% mkdir(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/SL_each/',expdir,exp,timeUnit,froidir));

for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~ ,tn ,subjn]=size(gdata);
    
    b=nan([length(rnames) length(lags) subjn]);
    r2=nan([length(rnames)  subjn]);
    % listener 10 in pieman did not go through the whole scanning (only 287TRs were obtained)
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            
            roi_voxn=size(gdata,1);
            keptT_s= max(-min(lags),0)+crop_start+1;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            for si=1:subjn;
                y=gdata(:,:,si);
                y=y(:,keptT);
                y=y(:);
                
                for li=1:length(lags);
                    X(:,:,li)=data(:,keptT+lags(li));
                end
                
                X=reshape(X,roi_voxn*length(keptT),length(lags));
                % centralize X
                X=X-mean(X);
                
                % include intercept
                [btemp,~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                b(ri,:,si)=btemp(2:end);
                r2(ri,si)=stats(1);
                clear X
            end
        end
    end
    
    save(sprintf('%s/%s/fmri/pattern/regression/%s/roi/%s/SL_each/lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'r2','b','rnames');
    
    clear b
end


