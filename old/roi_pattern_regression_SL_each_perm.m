function roi_tr_pattern_regression_SLeach_perm(perm)

disp(perm)
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=1:4;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
        [~ ,tn ,subjn]=size(gdata);
        
        b=nan([length(rnames) length(lags)+1 subjn]);
        % listener 10 in pieman did not go through the whole scanning (only 287TRs were obtained)
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
                
                roi_voxn=size(gdata,1);
                
                keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
                keptT_e=min(tn,tn-max(lags));
                keptT=keptT_s:keptT_e;
                
                if perm== 0; % speaker as predictor
                    data_perm=data;
                    gdata_perm=gdata;
                else
                    gdata_perm=gdata;
                    data_perm=gdata(:,:,perm);
                    gdata_perm(:,:,perm)=NaN;
                end
                
                for si=1:subjn;
                    if si ~= perm;
                        y=gdata_perm(:,:,si);
                        y=y(:,keptT);
                        y=y(:);
                        
                        for li=1:length(lags);
                            X(:,:,li)=data_perm(:,keptT+lags(li));
                        end
                        
                        X=reshape(X,roi_voxn*length(keptT),length(lags));
                        % centralize X
                        X=X-mean(X);
                        
                        % include intercept
                        [b(ri,:,si),~,~,~,~]=regress(y,[ones(size(X,1),1) X]);
                        
                        clear X
                    end
                end
            end
        end
        
        save(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/SLeach/perm/regression_SLeach_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b');
        
        clear b
    end
end

