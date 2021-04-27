%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

lags_tested={-40:40,-10:10,-15:15, -20:20};
permN=10000;

for ei=2:10;%
    exp=exp_parameters.experiments{ei};
    
    for lagi=4%:length(lags_tested);
        lags=lags_tested{lagi};
        
        rz=[];
        for sdi=1:length(networks);
            seed=networks{sdi};
            f=ls([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/' seed '.mat']);
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/' f ],'networks','r','keptT');
            rz(sdi,:,:,:)=atanh(r);
        end
        rzm=nanmean(rz,4);
        [~,~,tn]=size(rzm);
        tmid_real=(tn-1)/2+1;
        
        tmid_shifts=1:tn;
        tmid_shifts=tmid_shifts((tmid_shifts+min(lags))>=1 & (tmid_shifts+max(lags))<=tn);
        peaks_shift=[];
        
        for perm=1:permN;
            ti=randi(length(tmid_shifts));
            tmid_shift=tmid_shifts(ti);
            [peaks_shift(:,:,perm),lagi]=max(rzm(:,:,tmid_shift+lags),[],3);
        end

        rzm_temp=rzm(:,:,tmid_real+lags);
          rzm_z=zscore(rzm,0,3);
        rzm_z_temp=rzm_z(:,:,tmid_real+lags);
        
        p=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>rzm_temp,4);
    
        
        peaks=nan(size(rzm,1),size(rzm,2));
        peakLags=nan(size(rzm,1),size(rzm,2));
        p_peaks=nan(size(rzm,1),size(rzm,2));;
        peaks_z=nan(size(rzm,1),size(rzm,2));
        
        for sdi=1:length(networks);
            for tgi=1:length(networks);
                temp=squeeze(rzm_temp(sdi,tgi,:));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                
                if ~isempty(pks);
                    pk=max(pks);
                    [lagi]=find(temp==pk);
                    peakLags(sdi,tgi)=lags(lagi);
                    peaks(sdi,tgi)=pk;
                    p_peaks(sdi,tgi)=p(sdi,tgi,lagi);
                    peaks_z(sdi,tgi)=rzm_z_temp(sdi,tgi,lagi);
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaks' ],...
            'networks','rzm','rz','lags','keptT','p','peaks_z','peakLags','peaks','p_peaks');
    end
end

