t%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-40:40,-10:10,-15:15, -20:20};
permN=10000;

for ei=[4]
    exp=experiments{ei};
    
    for lagi=4%:length(lags_tested);
        lags=lags_tested{lagi};
        
        bs=[];
        for tgi=1:length(networks);
           target=networks{tgi};
            f=ls([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/LL_leave1out/target_' target '.mat']);
            load([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/LL_leave1out/' f ],'networks','b','keptT');
            bs(:,tgi,:,:)=b;
            
        end
        bm=nanmean(bs,4);
        [~,~,tn]=size(bm);
        t_real=(tn-1)/2+1;
        
        ts_shift=1:tn;
        ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
        peaks_shift=[];
        
        for perm=1:permN;
            ti=randi(length(ts_shift));
            t_shift=ts_shift(ti);
            [peaks_shift(:,:,perm),lagi]=max(bm(:,:,t_shift+lags),[],3);
        end
        
        bm=bm(:,:,t_real+lags);
        p=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>bm,4);
        pfwe=p*length(p(:))/length(lags);
        
        peaks=nan(size(bm,1),size(bm,2));
        peakLags=nan(size(bm,1),size(bm,2));
        peaks_pfwe=nan(size(bm,1),size(bm,2));
        peakLags_pfwe=nan(size(bm,1),size(bm,2));
        peaks_p=nan(size(bm,1),size(bm,2));
        peakLags_p=nan(size(bm,1),size(bm,2));
        
        for sdi=1:length(networks);
            for ni=1:length(networks);
                temp=squeeze(bm(sdi,ni,:));
                [pks, locs]=findpeaks(temp);
                locs=locs(pks>0);
                pks=pks(pks>0);
                if ~isempty(pks);
                    [~,loci]=max(pks);
                    peakLags(sdi,ni)=lags(locs(loci));
                    peaks(sdi,ni)=temp(locs(loci));
                    
                    if p(sdi,ni,locs(loci))<0.05;
                        peakLags_p(sdi,ni)=lags(locs(loci));
                        peaks_p(sdi,ni)=temp(locs(loci));
                    end
                    
                    if pfwe(sdi,ni,locs(loci))<0.05;
                        peakLags_pfwe(sdi,ni)=lags(locs(loci));
                        peaks_pfwe(sdi,ni)=temp(locs(loci));
                    end
                end
                
                
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
            'networks','bm','bs','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe');
    end
end

