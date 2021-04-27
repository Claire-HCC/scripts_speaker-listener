%% find the peak nearest to lag 0 instead of the absolute peak
clear all

% loc='cluster';
set_parameters;

timeUnit='tr' ;
froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
lags=-15:15;

rzms=[];
rzms_timeReversed=[];

for ei=1:8;
    exp=exp_parameters.experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/networks2networks_lag-15-15_timeReversed_peaks' ],'rzm','rzm_timeReversed','keptT','networks');
    
    [~,~,tn]=size(rzm);
    tmid=(tn-1)/2+1;
    
    rzms(:,:,:,ei)=rzm(:,:,tmid+[-100:100]);
    rzms_timeReversed(:,:,:,ei)=rzm_timeReversed(:,:,tmid+[-100:100]);
end

rzmsm=nanmean(rzms,4);
rzmsm_timeReversed=nanmean(rzms_timeReversed,4);

% paired t-test
tmid=(size(rzms,3)-1)/2+1;
p_time=[];
z_time=[];
% paired t-test
for sdi=1:length(networks);
    for tgi=1:length(networks);
        null=(squeeze(rzmsm_timeReversed(sdi,tgi,:)));
        null_m=mean(null);
        null_std=std(null);
        for lagi=1:length(lags);
            r_real=squeeze(rzmsm(sdi,tgi,tmid+lags(lagi)));
            [~,p_time(sdi,tgi,lagi),~,z_time(sdi,tgi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
        end
    end
end
[~,~,pfdr_time]=fdr(p_time(:),.05);
pfdr_time=reshape(pfdr_time,size(p_time));

p_exps=[];
t_exps=[];
peakLags=[];
for sdi=1:length(networks);
    for tgi=1:length(networks);
        for lagi=1:length(lags);
            r_real=squeeze(rzms(sdi,tgi,tmid+lags(lagi),:))';
            [~,p_exps(sdi,tgi,lagi),~,stats] = ttest(r_real,0,'tail','right');
            t_exps(sdi,tgi,lagi)=stats.tstat;
        end
    end
end

[~,~,pfdr_exps]=fdr(p_exps(:),.05);
pfdr_exps=reshape(pfdr_exps,size(p_exps));

peaks=nan(size(rzm,1),size(rzm,2));
peakLags=nan(size(rzm,1),size(rzm,2));
p_time_peaks=nan(size(rzm,1),size(rzm,2));
pfdr_time_peaks=nan(size(rzm,1),size(rzm,2));
p_exps_peaks=nan(size(rzm,1),size(rzm,2));
pfdr_exps_peaks=nan(size(rzm,1),size(rzm,2));
t_peaks=nan(size(rzm,1),size(rzm,2));

npeaks=nan(size(rzm,1),size(rzm,2));
npeakLags=nan(size(rzm,1),size(rzm,2));
% p_time_npeaks=nan(size(rzm,1),size(rzm,2));
% pfdr_time_npeaks=nan(size(rzm,1),size(rzm,2));
% p_exps_npeaks=nan(size(rzm,1),size(rzm,2));
% pfdr_exps_npeaks=nan(size(rzm,1),size(rzm,2));
% t_npeaks=nan(size(rzm,1),size(rzm,2));

for sdi=1:length(networks);
    for tgi=1:length(networks);
        
        temp=squeeze(rzmsm(sdi,tgi,tmid+lags));
        [pks]=findpeaks(temp);
        pks=pks(pks>0);
        
        if ~isempty(pks);
            pk=max(pks);
            [lagi]=find(temp==pk);
            peakLags(sdi,tgi)=lags(lagi);
            peaks(sdi,tgi)=pk;
            p_time_peaks(sdi,tgi)=p_time(sdi,tgi,lagi);
            pfdr_time_peaks(sdi,tgi)=pfdr_time(sdi,tgi,lagi);
            p_exps_peaks(sdi,tgi)=p_exps(sdi,tgi,lagi);
            pfdr_exps_peaks(sdi,tgi)=pfdr_exps(sdi,tgi,lagi);
            %    t_peaks(sdi,tgi)=t(sdi,tgi,lagi);
        end
        
        temp=-squeeze(rzmsm(sdi,tgi,tmid+lags));
        [pks]=findpeaks(temp);
        pks=pks(pks>0);
        if ~isempty(pks);
            pk=-max(pks);
            [lagi]=find(temp==max(pks));
            npeakLags(sdi,tgi)=lags(lagi);
            npeaks(sdi,tgi)=pk;
%             p_time_npeaks(sdi,tgi)=p_time(sdi,tgi,lagi);
%             pfdr_time_npeaks(sdi,tgi)=pfdr_time(sdi,tgi,lagi);
%             p_exps_npeaks(sdi,tgi)=p_exps(sdi,tgi,lagi);
%             pfdr_exps_npeaks(sdi,tgi)=pfdr_exps(sdi,tgi,lagi);
            %      t_npeaks(sdi,tgi)=t(sdi,tgi,lagi);
        end
    end
end

exp='crossExps';
mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/' ]);
rzm=rzmsm;
pfdr_peaks=pfdr_exps_peaks;
pfdr_peaks(:,:,2)=pfdr_time_peaks;
pfdr_peaks=max(pfdr_peaks,[],3);
% pfdr_npeaks=pfdr_exps_npeaks;
% pfdr_npeaks(:,:,2)=pfdr_time_npeaks;
% pfdr_npeaks=max(pfdr_npeaks,[],3);

save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
    'networks','rzms','rzms_timeReversed','rzmsm','rzmsm_timeReversed','lags','keptT','peakLags','peaks','rzm','pfdr_peaks','pfdr_exps_peaks','pfdr_time_peaks','tmid','npeaks');



