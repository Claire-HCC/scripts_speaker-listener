
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
[~,rnames,~] = cellfun(@fileparts,cellstr(ls([expdir '/roi_mask/'  froidir '/mat/*.mat'])),'UniformOutput',0);

tic % 15 min

lags=-40:40;
for ei=2%:2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                for s=1:size(gdata,3);
                    
                    self=gdata(:,:,s);
                    othersi=1:size(gdata,3);
                    othersi(othersi==s)=[];
                    others=nanmean(gdata(:,:,othersi),3);
                    
                    sim_SL(:,:,s)=lagcorr_claire(self,data,lags);
                    sim_LL(:,:,s)=lagcorr_claire(self,others,lags);
                end
                
                sim=0.5*log((1+sim_SL)./(1-sim_SL));
                save([expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/simz_SL_' rname ],'sim','lags');
                sim=0.5*log((1+sim_LL)./(1-sim_LL));
                save([expdir '/' exp '/fmri/sim/' timeUnit '/roi/' froidir '/simz_LL_' rname ],'sim','lags');
            end
        end
    end
end

toc
