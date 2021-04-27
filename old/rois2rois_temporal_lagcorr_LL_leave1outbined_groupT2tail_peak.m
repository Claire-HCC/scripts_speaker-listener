%% ttest across subject
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};
binSize=250;

for ei=3:4;%1:4;%1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' rnames{1} '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        [~,tn,~,listenerN]=size(r);
        bis=find(~isnan(squeeze(r(1,:,1,1))));
        
        for bi = [min(bis) max(bis)];
            
            rz=nan([length(rnames) length(rnames) length(lags) listenerN]);
            for sdi=1:length(rnames);
                seed=rnames{sdi};
                if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                    load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
                    rz(sdi,:,:,:)=atanh(squeeze(r(:,bi,:,:)));
                end
            end
            
            rzm=nanmean(rz,4);
            p=nan(size(rzm));
            t=p;
            for sdi=1:length(rnames);
                for ri=1:length(rnames);
                    [~,p(sdi,ri,:),~,stats]=ttest(squeeze(rz(sdi,ri,:,:))',0,'tail','both');
                    t(sdi,ri,:)=stats.tstat;
                end
            end
            
            pfwe=p*(sum(~isnan(p(:))));
            
            pfdr=nan(size(p(:)));
            [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
            pfdr=reshape(pfdr,size(p));
            
            peaks=nan([length(rnames) length(rnames) ]);
            peakLags=peaks;
            peaks_pfwe=nan([length(rnames) length(rnames) ]);
            peakLags_pfwe=peaks_pfwe;
            peaks_pfdr=nan([length(rnames) length(rnames) ]);
            peakLags_pfdr=peaks_pfdr;
            
            rz_tempfwe=rzm;
            rz_tempfwe(pfwe>.05)=NaN;
            
            rz_tempfdr=rzm;
            rz_tempfdr(pfdr>.05)=NaN;
            
            for sdi=1:length(rnames);
                for ri=1:length(rnames);
                    if sum(isnan(rzm(sdi,ri,:))==0);
                        
                        [~, peakLagi]=max(abs(rzm(sdi,ri,:)),[],3);%max(abs(rzm(ri,:)),[],2);
                        peakLags(sdi,ri)=(lags(peakLagi));
                        peaks(sdi,ri)=rzm(sdi,ri,peakLagi);
                        
                        if min(pfdr(sdi,ri,:))<.05;
                            [~, peakLagi]=max(abs(rz_tempfdr(sdi,ri,:)),[],3);;%max(abs(rz_tempfdr(sdi,ri,:)),[],2);
                            peakLags_pfdr(sdi,ri)=(lags(peakLagi));
                            peaks_pfdr(sdi,ri)=rz_tempfdr(sdi,ri,peakLagi);
                        end
                        
                        if min(pfwe(sdi,ri,:))<.05;
                            [~, peakLagi]=max(abs(rz_tempfwe(sdi,ri,:)),[],3);%max(abs(rz_tempfwe(sdi,ri,:)),[],2);
                            peakLags_pfwe(sdi,ri)=(lags(peakLagi));
                            peaks_pfwe(sdi,ri)=rz_tempfwe(sdi,ri,peakLagi);
                        end
                    end
                end
            end
            
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/rois2rois_binSize' num2str(binSize) '_bin' num2str(bi) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
               'rnames','t','r','rzm','lags','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
        end
    end
end
