clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
win=-5:15;
types={'wd','sn','pr'};
crop_start=25;
crop_end=20;

xBF.dt=1.5;
xBF.name='hrf (with time and dispersion derivatives)';
bf = spm_get_bf(xBF)

for ei=1;% [4 11 9 10];
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    mkdir([expdir '/' exp '/fmri/erp/network/' froidir ]);
    load([expdir exp '/sound/onsets.mat'],'onsets_wd','onsets_sn','onsets_pr','onsetsVector_wd','onsetsVector_sn','onsetsVector_pr');
    onsets_wd(ismember(onsets_wd,onsets_sn) )=[];
    onsets_wd(ismember(onsets_wd,onsets_pr) )=[];
    onsets_sn(ismember(onsets_sn,onsets_pr) )=[];
    
    load('Y:\claire\speaker-listener\sherlock\sound\sherlock_listener_audenv.mat','aud')
    aud=zscore(aud);
    
    for typei=1:3;
        
        type=types{typei};
        eval(['onsets=onsets_' type ';'])
        erps=nan(length(networks),length(win),length(onsets));
        erps_aud= nan(1,length(win),length(onsets));
        erps_onsetsWD= nan(1,length(win),length(onsets));
        
        for sdi=1:length(networks);
                        seed=networks{sdi};
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
                        % keep only voxel with peakLag 0
            
%                         [~,tn,listenerN]=size(gdata);
%                         keptT=(crop_start+1):(tn-crop_end);
%             
%                         gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
%                         gdata(:,:,subjects_excluded{ei})=NaN;
%                         gdata=nanmean(gdata,1);
%                         gdata(:,keptT,:)=zscore(gdata(:,keptT,:),0,2);
%                         gdata=nanmean(gdata,3);
            gdata=nan(1,size(gdata,2));
            load([expdir '/' exp '/fmri/erp/network/' froidir '/regress_text.mat'   ],'resid','keptT');
            gdata(:,keptT,:)=nanmean(zscore(resid(sdi,:,:),0,2),3);
            
            for oi=1:length(onsets);
                onset=round(onsets(oi));
                trs=onset+win;
                if sum(~ismember(trs,keptT))==0 & (onset+max(win))<=max(keptT);
                    
                    erps(sdi,:,oi)=gdata(trs);
                    
                    erps_aud(1,:,oi)=aud(trs);
                    erps_onsetsWD(1,:,oi)= onsetsVector_wd(trs);
                end
            end
        end
        
        
        save([expdir '/' exp '/fmri/erp/network/' froidir '/erp_' type '.mat'   ],'erps','networks','keptT','win','types','erps_aud','erps_onsetsWD');
    end
end



