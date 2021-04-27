%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='fc_cluster';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% networks=unique(table2array(roi_table(:,2)));
lags_tested={-40:40,-10:10,-15:15, -20:20};
permN=10000;
networks={'network1','network2','network3','network4','network5','network6'};

for ei=[7 8];%[1 13];%[1 2 4 11 12  13]
    exp=exp_parameters.experiments{ei};
    
    for lagi=4%:length(lags_tested);
        lags=lags_tested{lagi};
        
        rz=[];
        for sdi=1:length(networks);
            seed=networks{sdi};
            f=ls([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_gg/' seed '.mat']);
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_gg/' f ],'networks','r','keptT');
            rz(sdi,:,:,:)=atanh(r);
            
        end
        rzm=nanmean(rz,4);
        [~,~,tn]=size(rzm);
        t_real=(tn-1)/2+1;
        
        ts_shift=1:tn;
        ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
        peaks_shift=[];
        
        for perm=1:permN;
            ti=randi(length(ts_shift));
            t_shift=ts_shift(ti);
            [peaks_shift(:,:,perm),lagi]=max(rzm(:,:,t_shift+lags),[],3);
        end
        
        rz=rz(:,:,t_real+lags,:);
        rzm=rzm(:,:,t_real+lags);
        p=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>rzm,4);
        pfwe=p*length(p(:))/length(lags);
        
        peaks=nan(size(rzm,1),size(rzm,2));
        peakLags=nan(size(rzm,1),size(rzm,2));
        peaks_pfwe=nan(size(rzm,1),size(rzm,2));
        peakLags_pfwe=nan(size(rzm,1),size(rzm,2));
        peaks_p=nan(size(rzm,1),size(rzm,2));
        peakLags_p=nan(size(rzm,1),size(rzm,2));
        
        for sdi=1:length(networks);
            for ni=1:length(networks);
                temp=squeeze(rzm(sdi,ni,:));
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
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_gg/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
            'networks','rzm','rz','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe');
    end
end

