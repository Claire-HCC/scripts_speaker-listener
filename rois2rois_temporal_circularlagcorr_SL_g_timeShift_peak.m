%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,-15:15. -20:20};
permN=10000;

for ei=2;%
    exp=experiments{ei};
    
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};

        for sdi=1:length(rnames);
            seed=rnames{sdi};
            f=([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/' seed '.mat']);
            
            if exist(f);
                load(f,'r','keptT');
                if sdi==1;
                    rz=nan(size(r,1),size(r,1),size(r,2));
                end
                rz(sdi,:,:)=atanh(r);
            end
        end
        rzm=rz;
        [~,~,tn]=size(rzm);
        t_real=(tn-1)/2+1;
        
        ts_shift=1:tn;
        ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
        peaks_shift=nan(length(rnames),length(rnames),length(lags));
        
        for perm=1:permN;
            ti=randi(length(ts_shift));
            t_shift=ts_shift(ti);
            [peaks_shift(:,:,perm),~]=max(rzm(:,:,t_shift+lags),[],3);
        end
        
        rz=rz(:,:,t_real+lags,:);
        rzm=rzm(:,:,t_real+lags);
        p=[];
        for lagi=1:length(lags);
               p(:,:,lagi)=mean(peaks_shift>squeeze(rzm(:,:,lagi)),3);
        end

       pfwe=p*length(p(:))/length(lags);
      
        peaks=nan(size(rzm,1),size(rzm,2));
        peakLags=nan(size(rzm,1),size(rzm,2));
        peaks_pfwe=nan(size(rzm,1),size(rzm,2));
        peakLags_pfwe=nan(size(rzm,1),size(rzm,2));
        peaks_p=nan(size(rzm,1),size(rzm,2));
        peakLags_p=nan(size(rzm,1),size(rzm,2));
        
        for sdi=1:length(rnames);
            for ni=1:length(rnames);
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
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
            'rnames','rzm','rz','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','peaks_p','peakLags_p');
        
    end
end
