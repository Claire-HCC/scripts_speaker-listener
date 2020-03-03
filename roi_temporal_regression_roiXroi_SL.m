
% loc='cluster';
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=[1 2 4];%1:4;%
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_roisM' ],'data');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_roisM' ],'gdata');
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        [~ ,tn ,subjn]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        b=nan([length(rnames) length(rnames) length(lags)+1 subjn]);
        r2=nan([length(rnames) length(rnames) subjn ]);
        F=nan([length(rnames) length(rnames) subjn]);
        p=nan([length(rnames) length(rnames) subjn]);
        r=nan([length(rnames) length(rnames) tn subjn]);
        
        % listener 10 in pieman did not go through the whole scanning (only 287TRs were obtained)
        
        for ri_s=1:length(rnames);
            for ri_l=1:length(rnames);
                if nansum(gdata(ri_l,:,1))~=0 & nansum(data(ri_s,:))~=0
                    for si=1:subjn;
                        
                        y=gdata(ri_l,:,si);
                        y=y(:,keptT);
                        y=y(:);
                        
                        for li=1:length(lags);
                            X(:,:,li)=data(ri_s,keptT+lags(li));
                        end
                        
                        X=reshape(X,length(keptT),length(lags));
                        % centralize X
                        X=X-mean(X);
                        
                        % include intercept
                        [b(ri_s,ri_l,:,si),~,r(ri_s,ri_l,keptT),~,stats]=regress(y,[ones(size(X,1),1) X]);
                        
                        r2(ri_s,ri_l,si)=stats(1);
                        F(ri_s,ri_l,si)=stats(2);
                        p(ri_s,ri_l,si)=stats(3);
                        
                        clear X
                    end
                end
            end
        end
        
        save(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLeach/regression_SLeach_roiXroi_lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'b','r2','r','p','F');
        clear b
    end
end


