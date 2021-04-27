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

for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'rnames','r','lags','keptT');
        rz_perm=atanh(r);
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        rz=atanh(r);
        
        p=min(mean(rz_perm>rz,3) , mean(rz_perm<rz,3));
        p(isnan(rz))=NaN;
        pfwe=p*(sum(~isnan(p(:))));
        
        pfdr=nan(size(p(:)));
        [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
        pfdr=reshape(pfdr,size(p));
        
        ris=find(sum(~isnan(rz),2)~=0);
        peak=nan([length(rnames) 1 ]);
        peakLags=peak;
        peaks_pfwe=nan([length(rnames) 1 ]);
        peakLags_pfwe=peaks_pfwe;
        peaks_pfdr=nan([length(rnames) 1 ]);
        peakLags_pfdr=peaks_pfdr;
        
        rz_tempfwe=rz;
        rz_tempfwe(pfwe>.025)=NaN;
        rz_tempfdr=rz;
        rz_tempfdr(pfdr>.025)=NaN;
        for i=1:length(ris);
            ri=ris(i);
            [~, peakLagi]=max(abs(rz(ri,:)),[],2);
            peakLags(ri,1)=(lags(peakLagi));
            peaks(ri,1)=rz(ri,peakLagi);
            
            if min(pfdr(ri,:))<.025;
                [~, peakLagi]=max(abs(rz_tempfdr(ri,:)),[],2);
                peakLags_pfdr(ri,1)=(lags(peakLagi));
                peaks_pfdr(ri,1)=rz_tempfdr(ri,peakLagi);
            end
            
            if min(pfwe(ri,:))<.025;
                [~, peakLagi]=max(abs(rz_tempfwe(ri,:)),[],2);
                peakLags_pfwe(ri,1)=(lags(peakLagi));
                peaks_pfwe(ri,1)=rz_tempfwe(ri,peakLagi);
            end
        end
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfwe) & peaks_pfwe>0),   peakLags_pfwe(~isnan( peaks_pfwe) & peaks_pfwe>0)]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_posPeakLags_pfwe.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfwe) & peaks_pfwe<0),   peakLags_pfwe(~isnan( peaks_pfwe) & peaks_pfwe<0)]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_negPeakLags_pfwe.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr) & peaks_pfdr>0),   peakLags_pfdr(~isnan( peaks_pfdr) & peaks_pfdr>0)]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_posPeakLags_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr) & peaks_pfdr<0),   peakLags_pfdr(~isnan( peaks_pfdr) & peaks_pfdr<0)]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_negPeakLags_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peaks_pfdr)),   peaks_pfdr(~isnan( peaks_pfdr))]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peaks_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peaks_pfwe)),   peaks_pfwe(~isnan( peaks_pfwe))]);
        nii.img(1,1,[1 2])=[-1 1];
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peaks_pfwe.nii']);
        
        save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],...
            'rnames','rz','r','lags','keptT','p','peaks','peakLags','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr');
    end
end


