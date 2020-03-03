
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

for ei=1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        %    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL_stats.mat' ],'sig_fdr');
        %   regression_mask=sig_fdr;
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=mean(r_perm>r,3);
        %  p(regression_mask==0,:)=NaN;
        p(isnan(r))=NaN;
        
        pfdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            pfdr=reshape(pfdr,size(p));
        end
        
        ris=find(sum(~isnan(r),2)~=0);
        peak=nan([length(rnames) 1 ]);
        peakLags=peak;
        [peak(ris) peakLagi]=max(r(ris,:),[],2);
        peakLags(ris,1)=(lags(peakLagi))';
        
        peak_pfdr=nan([length(rnames) 1 ]);
        peakLags_pfdr=peak_pfdr;
        r_temp=r;
        r_temp(pfdr>.05)=NaN;
        ris=find(sum(~isnan(r_temp),2)~=0);
        [peak_pfdr(ris) peakLagi_pfdr]=max(r_temp(ris,:),[],2);
        peakLags_pfdr(ris,1)=(lags(peakLagi_pfdr))';
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peak_pfdr)),   peak_pfdr(~isnan( peak_pfdr))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL_peak_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags)),   peakLags(~isnan(peakLags))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL_peakLags.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr)),   peakLags_pfdr(~isnan( peakLags_pfdr))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/'  froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL_peakLags_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL_peaks' ],'rnames','r','lags','keptT','p','pfdr','peak_pfdr','peak','peakLags','peakLags_pfdr');
    end
end

