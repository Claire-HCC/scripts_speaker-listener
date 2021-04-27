%% find the peak nearest to lag 0 instead of the absolute peak
clear all

% loc='cluster';
set_parameters;

timeUnit='tr' ;
froidir='restFc_isc30PercMasked_cluster6';
lags_tested={-40:40,-10:10,-15:15, -20:20};

networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

for ei=1%:11;
    exp=exp_parameters.experiments{ei};
    
    for lagi=4%:length(lags_tested);
        lags=lags_tested{lagi};
        
        rz=[];
        rz_perm=[];
        for sdi=1%:length(networks);
            seed=networks{sdi};
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/' seed ],'r','keptT');
            rz(sdi,:,:,:)=atanh(r);
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/perm/' seed '_permPhase' ],'r','keptT');
            rz_perm(sdi,:,:,:,:)=atanh(r);
        end
        
        [~,~,tn,listenerN]=size(rz);
        tmid=(tn-1)/2+1;
        rzm=nanmean(rz,4);
        rzm_perm=nanmean(rz_perm,4);
        p=nanmean(rzm_perm(:,:,tmid+lags,:)>rzm(:,:,tmid+lags),4);
        [~,~,pfdr]=fdr(p(:));
        pfdr=reshape(pfdr,size(p));
        
        peaks=nan(size(p,1),size(p,2));
        peakLags=nan(size(p,1),size(p,2));
        p_peaks=nan(size(p,1),size(p,2));
        pfdr_peaks=nan(size(p,1),size(p,2));
        
        
        for sdi=1:size(p,1)
            for ni=1:size(p,2)
                temp=squeeze(rzm(sdi,ni,tmid+lags));
                [pks, locs]=findpeaks(temp);
                locs=locs(pks>0);
                pks=pks(pks>0);
                
                if ~isempty(pks)
                    [~,loci]=max(pks);
                    peakLags(sdi,ni)=lags(locs(loci));
                    peaks(sdi,ni)=temp(locs(loci));
                    
                    p_peaks(sdi,ni)=p(sdi,ni,locs(loci));
                    pfdr_peaks(sdi,ni)=pfdr(sdi,ni,locs(loci));
                    
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_peaks' ],...
            'networks','rzm','rz','lags','keptT','p','pfdr','peakLags','peaks','p_peaks','pfdr_peaks');
        
    end
end

