
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags_tested={-4:4, -10:10};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','keptvox','p');
        pfdr=nan(size(p));
        [~,~,pfdr]=fdr(p);
        
        mat=nan(voxn,1);
        mat(keptvox(pfdr<.05))=r2(pfdr<.05);
        nii=mat2nii(mat);
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','keptvox','p','pfdr');
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_ftest_pfdr.nii']);
        

        
     
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
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfwe) & peaks_pfwe>0),   peakLags_pfwe(~isnan( peaks_pfwe) & peaks_pfwe>0)+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_posPeakLags_pfwe.nii']);
    end
end


