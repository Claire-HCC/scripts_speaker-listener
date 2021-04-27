%% find the peak nearest to lag 0 instead of the absolute peak
clear all

loc='mypc';
set_parameters;

timeUnit='tr' ;
froidir='restFc_isc30PercMasked_75Overlap_cluster6_audPausesResid';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

lags_tested={-40:40,-10:10,-15:15, -20:20};

for ei=2;%[1 2 6];%1:11;
    exp=exp_parameters.experiments{ei};
    
    for lagi=3%:length(lags_tested);
        lags=lags_tested{lagi};
        
        rz=[];
        rz_timeReversed=[];
        for sdi=1:length(networks);
            seed=networks{sdi};
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/' seed ],'r','keptT');
            rz(sdi,:,:,:)=atanh(r);
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/' seed '_timeReversed' ],'r','keptT');
            rz_timeReversed(sdi,:,:,:)=atanh(r);
        end
        rzm=nanmean(rz,4);
        rzm_timeReversed=nanmean(rz_timeReversed,4);
        
        [~,~,tn,listenerN]=size(rz);
        tmid=(tn-1)/2+1;
        p=[];
        z=[];
        % paired t-test
        for sdi=1:length(networks);
            for tgi=1:length(networks);
                null=(squeeze(rzm_timeReversed(sdi,tgi,:)));
                null_m=mean(null);
                null_std=std(null);
                for lagi=1:length(lags);
                    r_real=squeeze(rzm(sdi,tgi,tmid+lags(lagi)));
                    [~,p(sdi,tgi,lagi),~,z(sdi,tgi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
                end
            end
        end
        [~,~,pfdr]=fdr(p(:),.05);
        pfdr=reshape(pfdr,size(p));
        
        peaks=nan(size(rzm,1),size(rzm,2));
        peakLags=nan(size(rzm,1),size(rzm,2));
        p_peaks=nan(size(rzm,1),size(rzm,2));
        pfdr_peaks=nan(size(rzm,1),size(rzm,2));
        z_peaks=nan(size(rzm,1),size(rzm,2));
        
        npeaks=nan(size(rzm,1),size(rzm,2));
        npeakLags=nan(size(rzm,1),size(rzm,2));
        %         p_npeaks=nan(size(rzm,1),size(rzm,2));
        %         pfdr_npeaks=nan(size(rzm,1),size(rzm,2));
        %         z_npeaks=nan(size(rzm,1),size(rzm,2));
        
        for sdi=1:length(networks);
            for tgi=1:length(networks);
                temp=squeeze(nanmean(rz(sdi,tgi,tmid+lags,:),4));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                
                if ~isempty(pks);
                    pk=max(pks);
                    [lagi]=find(temp==pk);
                    peakLags(sdi,tgi)=lags(lagi);
                    peaks(sdi,tgi)=pk;
                    p_peaks(sdi,tgi)=p(sdi,tgi,lagi);
                    pfdr_peaks(sdi,tgi)=pfdr(sdi,tgi,lagi);
                    z_peaks(sdi,tgi)=z(sdi,tgi,lagi);
                end
                
                temp=-squeeze(nanmean(rz(sdi,tgi,tmid+lags,:),4));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                if ~isempty(pks);
                    pk=-max(pks);
                    [lagi]=find(temp==max(pks));
                    npeakLags(sdi,tgi)=lags(lagi);
                    npeaks(sdi,tgi)=pk;
                    %                     p_npeaks(sdi,tgi)=p(sdi,tgi,lagi);
                    %                     pfdr_npeaks(sdi,tgi)=pfdr(sdi,tgi,lagi);
                    %                     z_npeaks(sdi,tgi)=z(sdi,tgi,lagi);
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
            'networks','rzm','rzm_timeReversed','rz','lags','keptT','p','z','peakLags','peaks','p_peaks','z_peaks','pfdr_peaks','pfdr','tmid','npeakLags','npeaks');
    end
end




