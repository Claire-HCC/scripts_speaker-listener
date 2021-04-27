%% find the peak nearest to lag 0 instead of the absolute peak
clear all
loc='mypc'
set_parameters;
timeUnit='tr' ;

lags_tested={-40:40,-10:10,-15:15, -20:20};

froidir='isc_peak';
rois={'precuneus','Angular_L','HG_L'};
perc=0.3;

for ei=1:9;%1:10;
    exp=exp_parameters.experiments{ei};
    
    for lagi=3%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for ri=1:length(rois);
            roi=rois{ri};
            f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,exp,'tr',perc*100);
            load(f,'mask');
            keptvox_mask=find(mask);
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/'  roi '.mat'   ],'r','keptvox','keptT');
            rz=atanh(r(ismember(keptvox,keptvox_mask),:,:));
            rzm=nanmean(rz,3);
            
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/'  froidir '/LL_selfself/' roi '_timeReversed.mat'   ],'r','keptT');
            rz_timeReversed=atanh(r(ismember(keptvox,keptvox_mask),:,:));
            rzm_timeReversed=nanmean(rz_timeReversed,3);
            keptvox=intersect(keptvox,keptvox_mask);
            
            [~,tn,listenerN]=size(rz);
            tmid=(tn-1)/2+1;
            p=[];
            z=[];
            % paired t-test
            
            for vi=1:length(keptvox);
                null=(squeeze(rzm_timeReversed(vi,:)));
                null_m=mean(null);
                null_std=std(null);
                for lagi=1:length(lags);
                    r_real=squeeze(rzm(vi,tmid+lags(lagi)));
                    [~,p(vi,lagi),~,z(vi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
                end
            end
            
            [~,~,pfdr]=fdr(p(:),.05);
            pfdr=reshape(pfdr,size(p));
            pfwe=p*length(p(:));
            
            peaks=nan(size(rzm,1),1);
            peakLags=nan(size(rzm,1),1);
            p_peaks=nan(size(rzm,1),1);
            pfdr_peaks=nan(size(rzm,1),1);
            z_peaks=nan(size(rzm,1),1);
            
            npeaks=nan(size(rzm,1),1);
            npeakLags=nan(size(rzm,1),1);
            %             p_npeaks=nan(size(rzm,1),1);
            %             pfdr_npeaks=nan(size(rzm,1),1);
            %             z_npeaks=nan(size(rzm,1),1);
            
            for tgi=1:size(rzm,1);
                temp=squeeze(nanmean(rz(tgi,tmid+lags,:),3));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                
                if ~isempty(pks);
                    pk=max(pks);
                    [lagi]=find(temp==pk);
                    peakLags(tgi)=lags(lagi);
                    peaks(tgi)=pk;
                    p_peaks(tgi)=p(tgi,lagi);
                    pfdr_peaks(tgi)=pfdr(tgi,lagi);
                    z_peaks(tgi)=z(tgi,lagi);
                end
                
                temp=-squeeze(nanmean(rz(tgi,tmid+lags,:),3));
                [pks]=findpeaks(temp);
                pks=pks(pks>0);
                
                if ~isempty(pks);
                    pk=-max(pks);
                    [lagi]=find(temp==max(pks));
                    npeakLags(tgi)=lags(lagi);
                    npeaks(tgi)=pk;
                    %                     p_npeaks(tgi)=p(tgi,lagi);
                    %                     pfdr_npeaks(tgi)=pfdr(tgi,lagi);
                    %                     z_npeaks(tgi)=z(tgi,lagi);
                end
            end
            
            save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
                'keptvox','rzm','rzm_timeReversed','rz','lags','keptT','p_peaks','p','peakLags','peaks','z_peaks','pfdr','pfdr_peaks','npeakLags','npeaks','tmid');
            
                        % significant positive peak that is larger than npeak
                        sig=(pfdr_peaks<.05 & (abs(npeaks)<peaks | isnan(npeaks) ));
            
                        mat=nan(voxn,1);
                        mat(keptvox(sig))=peaks(sig);
                        nii=mat2nii(mat);
                        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fdr_peaks.nii' ]);
            
                        mat=nan(voxn,1);
                        mat(keptvox(sig))=peakLags
                        nii=mat2nii(mat);
                        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fdr_peakLags.nii' ]);
                        
                            mat=nan(voxn,1);
                        mat(keptvox)=peakLags;
                        nii=mat2nii(mat);
                        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peakLags.nii' ])

        end
    end
end


