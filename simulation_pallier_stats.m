loc='cluster';
set_parameters;
exp='simulation';

simN=30;
lags=-15:15;

filelist = dir(sprintf('%s/simulation/sim%03d/*_r.mat',expdir,simN));
filelist = filelist(~endsWith({filelist.name}, '_peaks.mat'));
filelist={filelist.name}';

for fi=1:length(filelist);
    if exist(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,[strrep(filelist{fi},'_r.mat','') '_peaks.mat']))==0;
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,filelist{fi}));
       
        rz=atanh(r-0.00000001);
        rz_timeReversed=atanh(r_timeReversed-0.00000001);
        rzm=nanmean(rz,4);
        rzm_timeReversed=nanmean(rz_timeReversed,4);
        
        [levn,~,tn,~]=size(rz);
        % paired t-test
        tmid=(tn-1)/2+1;
        p_time=[];
        z_time=[];
        % paired t-test
        for sdi=1:levn;
            for tgi=1:levn;
                null=(squeeze(rzm_timeReversed(sdi,tgi,:)));
                null_m=mean(null);
                null_std=std(null);
                for lagi=1:length(lags);
                    r_real=squeeze(rzm(sdi,tgi,tmid+lags(lagi)));
                    [~,p_time(sdi,tgi,lagi),~,z_time(sdi,tgi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
                end
            end
        end
        [~,~,pfdr_time]=fdr(p_time(:),.05);
        pfdr_time=reshape(pfdr_time,size(p_time));
        
        p_exps=[];
        t_exps=[];
        peakLags=[];
        for sdi=1:levn;
            for tgi=1:levn;
                for lagi=1:length(lags);
                    r_real=squeeze(rz(sdi,tgi,tmid+lags(lagi),:))';
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
        
        % it does not make sense to use the p-value for
        % nagative R with one-tail t-test
        npeaks=nan(size(rzm,1),size(rzm,2));
        npeakLags=nan(size(rzm,1),size(rzm,2));
        %                         p_time_npeaks=nan(size(rzm,1),size(rzm,2));
        %                         pfdr_time_npeaks=nan(size(rzm,1),size(rzm,2));
        %                         p_exps_npeaks=nan(size(rzm,1),size(rzm,2));
        %                         pfdr_exps_npeaks=nan(size(rzm,1),size(rzm,2));
        %                         t_npeaks=nan(size(rzm,1),size(rzm,2));
        
        for sdi=1:levn;
            for tgi=1:levn;
                
                temp=squeeze(rzm(sdi,tgi,tmid+lags));
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
                
                temp=-squeeze(rzm(sdi,tgi,tmid+lags));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                if ~isempty(pks) & sum(pks==max(pks))==1;
                    pk=-max(pks);
                    [lagi]=find(temp==max(pks));
                    npeakLags(sdi,tgi)=lags(lagi);
                    npeaks(sdi,tgi)=pk;
                    %                                     p_time_npeaks(sdi,tgi)=p_time(sdi,tgi,lagi);
                    %                                     pfdr_time_npeaks(sdi,tgi)=pfdr_time(sdi,tgi,lagi);
                    %                                     p_exps_npeaks(sdi,tgi)=p_exps(sdi,tgi,lagi);
                    %                                     pfdr_exps_npeaks(sdi,tgi)=pfdr_exps(sdi,tgi,lagi);
                    %      t_npeaks(sdi,tgi)=t(sdi,tgi,lagi);
                end
            end
        end
        
        pfdr_peaks=pfdr_exps_peaks;
        pfdr_peaks(:,:,2)=pfdr_time_peaks;
        pfdr_peaks=max(pfdr_peaks,[],3);
        %                         pfdr_npeaks=pfdr_exps_npeaks;
        %                         pfdr_npeaks(:,:,2)=pfdr_time_npeaks;
        %                         pfdr_npeaks=max(pfdr_npeaks,[],3);
        
        
        save(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,[strrep(filelist{fi},'_r.mat','') '_peaks.mat']),...
            'rzm','rzm_timeReversed','rz','lags','keptT','peakLags','peaks','npeakLags','npeaks','pfdr_peaks','pfdr_time_peaks','pfdr_exps_peaks','tmid');
        
        disp([strrep(filelist{fi},'_r.mat','') '_peaks.mat']);
        sig=(pfdr_time_peaks<.05 & pfdr_exps_peaks<.05 & (abs(npeaks)<peaks | isnan(npeaks) ));
        peakLags_pfdr=peakLags;
        peakLags_pfdr(sig==0)=NaN;
        peakLags_pfdr
    end
end
