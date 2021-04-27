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
        
        
        f=([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/r.mat']);
        load(f,'r','keptT');
        rz=atanh(r);
        
        [~,tn]=size(rz);
        t_real=(tn-1)/2+1;
        
        ts_shift=1:tn;
        ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
        peaks_shift=nan(length(rnames),length(lags));
        
        for perm=1:permN;
            ti=randi(length(ts_shift));
            t_shift=ts_shift(ti);
            [peaks_shift(:,perm),~]=max(rz(:,t_shift+lags),[],2);
        end
        
        rz=rz(:,t_real+lags,:);
        p=[];
        
        for lagi=1:length(lags);
            p(:,lagi)=mean(peaks_shift>squeeze(rz(:,lagi)),2);
        end
        
        pfwe=p*length(p(:))/length(lags);
        
        peaks=nan(size(rz,1),1);
        peakLags=nan(size(rz,1),1);
        peaks_pfwe=nan(size(rz,1),1);
        peakLags_pfwe=nan(size(rz,1),1);
        peaks_p=nan(size(rz,1),1);
        peakLags_p=nan(size(rz,1),1);
        
        
        for ri=1:length(rnames);
            temp=squeeze(rz(ri,:));
            [pks, locs]=findpeaks(temp);
            locs=locs(pks>0);
            pks=pks(pks>0);
            
            if ~isempty(pks);
                [~,loci]=max(pks);
                peakLags(ri,1)=lags(locs(loci));
                peaks(ri,1)=temp(locs(loci));
                
                if p(ri,locs(loci))<0.05;
                    peakLags_p(ri,1)=lags(locs(loci));
                    peaks_p(ri,1)=temp(locs(loci));
                end
                
                if pfwe(ri,locs(loci))<0.05;
                    peakLags_pfwe(ri,1)=lags(locs(loci));
                    peaks_pfwe(ri,1)=temp(locs(loci));
                end
            end
        end
        
        
        save([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/roi_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
            'rnames','rz','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','peaks_p','peakLags_p');
        
    end
end

