
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

for ei=[3 4 1 2];%1:4;%
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_roisM' ],'data');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_roisM' ],'gdata');
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        [~ ,tn ,subjn]=size(gdata);
        
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames) length(rnames) length(lags) subjn]);
        peakLag=nan([length(rnames) length(rnames) subjn]);
        
        % listener 10 in pieman did not go through the whole scanning (only 287TRs were obtained)
        
        for ri_s=1:length(rnames);
            for ri_l=1:length(rnames);
                if nansum(gdata(ri_l,:,1))~=0 & nansum(data(ri_s,:))~=0
                    for si=1:subjn;
                        
                        y=gdata(ri_l,keptT,si);
                        
                        x=data(ri_s,keptT);
                        [r(ri_s,ri_l,:,si) ,~]=lagcorr(y',x');
                        [~, lagi]=max(squeeze(r(ri_s,ri_l,:,si)));
                        peakLag(ri_s,ri_l,si)=lags(lagi);
                    end
                    clear x y
                end
            end
        end
    end
    mkdir(sprintf('%s/%s/fmri/temporal_isc/%s/roi/%s/SLeach/',expdir,exp,timeUnit,froidir));
    save(sprintf('%s/%s/fmri/temporal_isc/%s/roi/%s/SLeach/regression_SLeach_roiXroi_lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'r','peakLag');
    clear peakLag r
end



