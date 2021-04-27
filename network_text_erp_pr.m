clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

cols=jet(7);
cols=cols([1 3 4 5 6 7],:);

win=-40:40;
types={'wd','sn','pr'};
crop_start=25;
crop_end=20;

for ei=[1 2 6];
    exp=exp_parameters.experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/erp/network/' froidir ]);
    load([expdir exp '/sound/onsets.mat']);
    offsets_tr_wd(ismember(offsets_tr_wd,offsets_tr_sn) )=[];
    offsets_tr_wd(ismember(offsets_tr_wd,offsets_tr_pr) )=[];
    offsets_tr_sn(ismember(offsets_tr_sn,offsets_tr_pr) )=[];
    
    for typei=1:3;
        type=types{typei};
        
        if exist(['offsets_tr_' type ]);
            eval(['onsets=onsets_tr_' type ';'])
            eval(['offsets=offsets_tr_' type ';'])
            
            for sdi=1:length(networks);
                seed=networks{sdi};
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox')
                % keep only voxel with peakLag 0
                
                [~,tn,listenerN]=size(gdata);
                keptT=(crop_start+1):(tn-crop_end);
                
                gdata=nanmean(gdata,1);
                gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
                gdata(:,~ismember(1:tn,keptT),:)=NaN;
                gdata(:,keptT,:)=gdata(:,keptT,:)-nanmean(gdata(:,keptT,:),2);
                
                
                if sdi==1;
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
                        
                        
                    end
                end
            end
            save([expdir '/' exp '/fmri/erp/network/' froidir '/erp_' type  '.mat'   ],'erps_onset','erps_offset','networks','keptT','win','types');
        end
    end
end


