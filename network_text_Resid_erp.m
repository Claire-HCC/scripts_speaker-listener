clear all
close all
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6_audPrpauseResid';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

cols=jet(7);
cols=cols([1 3 4 5 6 7],:);

win=-40:40;
types={'wd','sn','pr'};
crop_start=25;
crop_end=20;

for ei=[6  ];
    exp=exp_parameters.experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/erp/network/' froidir ]);
    load([expdir exp '/sound/onsets.mat']);
    
    for typei=3;%1:3;
        type=types{typei};
        
        if exist(['offsets_tr_' type ]);
            eval(['onsets=onsets_tr_' type ';'])
            eval(['offsets=offsets_tr_' type ';'])
            
            for sdi=1:length(networks);
                seed=networks{sdi};
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox','Y','y')
                % keep only voxel with peakLag 0
                
                [~,tn,listenerN]=size(gdata);
                keptT=(crop_start+1):(tn-crop_end);
                
                gdata=nanmean(gdata,1);
                gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
                gdata(:,~ismember(1:tn,keptT),:)=NaN;
                gdata(:,keptT,:)=gdata(:,keptT,:)-nanmean(gdata(:,keptT,:),2);
                
 
                Y=nanmean(Y,1);
                Y(:,~ismember(1:tn,keptT),:)=NaN;
                Y(:,keptT,:)=Y(:,keptT,:)-nanmean(Y(:,keptT,:),2);
                
                y=nanmean(y,1);
                y(:,~ismember(1:tn,keptT),:)=NaN;
                y(:,keptT,:)=y(:,keptT,:)-nanmean(y(:,keptT,:),2);
                
                
                if sdi==1;
                    erpsy_onset=nan(length(networks),length(win), listenerN,length(onsets));
                    erpsy_offset=nan(length(networks),length(win), listenerN,length(offsets));
                    erpsY_onset=nan(length(networks),length(win), listenerN,length(onsets));
                    erpsY_offset=nan(length(networks),length(win), listenerN,length(offsets));
                    erps_onset=nan(length(networks),length(win), listenerN,length(onsets));
                    erps_offset=nan(length(networks),length(win), listenerN,length(offsets));
                end
                for oi=1:length(offsets);
                    offset=round(offsets(oi));
                    trs_offset=offset+win;
                    
                    onset=round(onsets(oi));
                    trs_onset=onset+win;
                    
                    if sum(~ismember(trs_onset,keptT))==0 & sum(~ismember(trs_offset,keptT))==0 & (offset+max(win))<=max(keptT) ;
                        
                        erps_onset(sdi,:,:,oi)=gdata(1,trs_onset,:);
                        erps_offset(sdi,:,:,oi)=gdata(1,trs_offset,:);
                        
                        erpsY_onset(sdi,:,:,oi)=Y(1,trs_onset,:);
                        erpsY_offset(sdi,:,:,oi)=Y(1,trs_offset,:);
                        
                        erpsy_onset(sdi,:,:,oi)=y(1,trs_onset,:);
                        erpsy_offset(sdi,:,:,oi)=y(1,trs_offset,:);
                        
                    end
                end
            end
            save([expdir '/' exp '/fmri/erp/network/' froidir '/erp_' type  '.mat'   ],'erps_onset','erps_offset','erpsY_onset','erpsY_offset','erpsy_onset','erpsy_offset','networks','keptT','win','types');
        end
    end
end


