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
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT','p');
        rz=atanh(r);
        rz(:,:,subjects_excluded{ei})=NaN;
        rzm=nanmean(rz,3);
        p=nan(size(rzm));
        t=p;
        for ri=1:length(rnames);
            [~,p(ri,:),~,stats]=ttest(squeeze(rz(ri,:,:))',0,'tail','right');
            t(ri,:)=stats.tstat;
        end
        
        pfwe=p*(sum(~isnan(p(:))));
        
        pfdr=nan(size(p(:)));
        [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
        pfdr=reshape(pfdr,size(p));
        
        ris=find(sum(~isnan(rzm),2)~=0);
        peak=nan([length(rnames) 1 ]);
        peakLags=peak;
        peaks_pfwe=nan([length(rnames) 1 ]);
        peakLags_pfwe=peaks_pfwe;
        peaks_pfdr=nan([length(rnames) 1 ]);
        peakLags_pfdr=peaks_pfdr;
        
        rz_tempfwe=rzm;
        rz_tempfwe(pfwe>.05)=NaN;
        rz_tempfdr=rzm;
        rz_tempfdr(pfdr>.05)=NaN;
        for i=1:length(ris);
            ri=ris(i);
            [~, peakLagi]=max(rzm(ri,:),[],2);;%max(abs(rzm(ri,:)),[],2);
            peakLags(ri,1)=(lags(peakLagi));
            peaks(ri,1)=rzm(ri,peakLagi);
            
            if min(pfdr(ri,:))<.05;
                [~, peakLagi]=max(rz_tempfdr(ri,:),[],2);%max(abs(rz_tempfdr(ri,:)),[],2);
                peakLags_pfdr(ri,1)=(lags(peakLagi));
                peaks_pfdr(ri,1)=rz_tempfdr(ri,peakLagi);
            end
            
            if min(pfwe(ri,:))<.05;
                [~, peakLagi]=max(rz_tempfwe(ri,:),[],2);% max(abs(rz_tempfwe(ri,:)),[],2);
                peakLags_pfwe(ri,1)=(lags(peakLagi));
                peaks_pfwe(ri,1)=rz_tempfwe(ri,peakLagi);
            end
        end
        
        % +0.0000000001, otherwise, xhview won't show 0 peak lag.
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfwe) & peaks_pfwe>0),   peakLags_pfwe(~isnan( peaks_pfwe) & peaks_pfwe>0)+0.0000000001]);
        % make sure that there are both positive and negative values in the
        % image. Otherwise, xjvoew won't use the right colorbar;
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags_pfwe.nii']);
        
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfwe) & peaks_pfwe<0),   peakLags_pfwe(~isnan( peaks_pfwe) & peaks_pfwe<0)+0.0000000001]);
        %         nii.img(1,1,[1 2])=[-1 1];
        %         save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_negPeakLags_pfwe.nii']);
        delete([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_negPeakLags_pfwe.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr) & peaks_pfdr>0),   peakLags_pfdr(~isnan( peaks_pfdr) & peaks_pfdr>0)+0.0000000001]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_posPeakLags_pfdr.nii']);
        
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr) & peaks_pfdr<0),   peakLags_pfdr(~isnan( peaks_pfdr) & peaks_pfdr<0)+0.0000000001]);
        %         nii.img(1,1,[1 2])=[-1 1];
        %         save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_negPeakLags_pfdr.nii']);
        delete([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_negPeakLags_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peaks_pfdr) ),   peaks_pfdr(~isnan( peaks_pfdr))]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_peaks_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peaks_pfwe) ),   peaks_pfwe(~isnan( peaks_pfwe))]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_groupT_peaks_pfwe.nii']);
        
        save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
            'rnames','rzm','r','lags','keptT','p','peaks','peakLags','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr');
    end
end


