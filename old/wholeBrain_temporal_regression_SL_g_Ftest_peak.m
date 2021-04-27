
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags_tested={-4:4, -10:10};

for ei=1;%1:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'pfdr','b');
      %  load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r');
        b=b(:,2:end);
        
        peaks_pfdr=nan([length(keptvox) 1 ]);
        peakLags_pfdr=peaks_pfdr;
        
        vis=find(pfdr<.05);
        rz_tempfdr=b;
        rz_tempfdr(pfdr>.05,:)=NaN;
        for vi=1:length(keptvox);
            if pfdr(vi)<.05;
                [~, peakLagi]=max(abs(rz_tempfdr(vi,:)),[],2);
                peakLags_pfdr(vi,1)=(lags(peakLagi));
                peaks_pfdr(vi,1)=rz_tempfdr(vi,peakLagi);
            end
            
        end
        
        mat=nan([voxn 1]);
        mat(keptvox(peaks_pfdr>0))=peakLags_pfdr(peaks_pfdr>0);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_posPeakLags_pfdr.nii']);
        
        mat=nan([voxn 1]);
        mat(keptvox(peaks_pfdr<0))=peakLags_pfdr(peaks_pfdr<0);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_negPeakLags_pfdr.nii']);
        
        mat=nan([voxn 1]);
        mat(keptvox)=peakss_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peaks_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],...
            'rnames','b','r','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr');
    end
end


