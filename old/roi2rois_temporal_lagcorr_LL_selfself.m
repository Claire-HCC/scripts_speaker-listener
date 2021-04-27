function roi2rois_temporal_lagcorr_LL_selfself(sdi)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
% seeds=rnames(~ismember(1:length(rnames),[1 14 11]));%{'HG_L','vPCUN','pANG_L'};
seeds=rnames;
seed=seeds{sdi};
crop_start=10;
lags_tested={-10:10, -40:40, -60:60};

for ei=3:4;
    exp=experiments{ei};
  %  rmdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_selfself/'],'s');
  %  mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_selfself/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};

        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=zscore(nanmean(gdata_seed,1),0,2);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames)  length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=zscore(nanmean(gdata,1),0,2);
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=gdata(:,keptT,si);
                    x=gdata_seed(:,keptT,si);
                    
                    r(ri,:,si)=lagcorr(x',y',lags);
                    
                end
            end
        end
        mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_selfself/']);
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_selfself/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear r
        % end
    end
end

