clear all
close all
loc='mypc'
set_parameters

role='listenerZscoreMean_g2';
relation='LL';
p_thr=0.01;
lags=[-10:10]; %scan

mni=load_nii([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);

rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
rois_selected={'HG_L','PMC_L','vPCUN','HG_R','vmPFC'};
ris=(find(contains(rtable.region,rois_selected)));

mask=  load_nii([expdir 'roi_mask/gray_matter_mask.nii']);
mask=mask.img(:);

for ei=1%:2;%1:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat//roi/' exp '_listenerZscoreMean_g1_rois.mat']);
    %   load([expdir experiments{ei} '/fmri/mat/roi/merlin_speaker_rois.mat']);
    speaker=data;
    
    for ri_selected=2:3;%2:length(ris);
        
        ri=ris(ri_selected);
        rname=rois_selected{ri_selected};
        
        speaker_roi=speaker(ri,:);
        speaker_roi=zscore(speaker_roi')';
        
        subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*' role '*.mat']));
        %  subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))))
        
        rid=rtable.id{ri};
        roi=load_nii([expdir '/roi_mask/nifti/' rid '.nii']);
        
        % compute listener-speaker isfc for each listener
        for si = 1%:size(listeners,3);
            % choose listener and take zscore of timeseries in each voxel.
            [filepath,subj,ext] = fileparts(subjects{si});
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj '.mat']);
            listener=data;
            
            frois=load_nii([expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags)) '_' rname '_peakT_p' num2str(p_thr) 'FDR_clusterMask.nii']);
            frois=frois.img;
            froi_labels=unique(frois(frois~=0));
            froiN=length(froi_labels);
            
            for fri=1:froiN;
                froi_label=froi_labels(fri);
                froi_mask=(frois==froi_label);
                listener_froi=mean(listener(froi_mask,:));
                
                mask_properties = bwconncomp(froi_mask);
                
                % Note though that x coordinate corresponds to column in MATLAB, so if you want to use the centroid to index you should reverse the first two coordinates
                mask_centroid=regionprops(mask_properties ,'Centroid');
                mask_centroid=mask_centroid.Centroid([2 1 3]);
                
                try
                    centroid_mni=cor2mni(mask_centroid);
                    [onelinestructure, cellarraystructure]= cuixuFindStructure(centroid_mni);
                    frname=strrep(strrep(cellarraystructure{6},'_',' '),' (aal)','');
                catch
                    frname='undefined';
                end
                
                %   lagcc_temp= lagcorr_claire(speaker_roi',listener_froi',lags)';
                lagcc_temp= lagcorr_claire(speaker_roi',listener_froi',[-40:40])';
                lag=sign(froi_label)*floor(abs(froi_label/100));
                
                fsize=[9 18];
                figure('papersize',fsize,'paperposition',[1 1 fsize],'position',[1 1 fsize],'unit','centimeter');
                
                subplot(3,2,1);
                overlay=mni.img(round(mask_centroid(1)),:,:);
                overlay((froi_mask(round(mask_centroid(1)),:,:)))=10000;
                imagesc(imrotate(squeeze(overlay),90));
                colormap('gray');
                title(frname);
                
                subplot(3,2,2);
                %     findpeaks(lagcc_temp,lags,'MinPeakHeight',0.15);
                findpeaks(lagcc_temp,[-40:40],'MinPeakHeight',0.15);
                hold on
                ylim([-0.4 0.4])
                plot([-40:40],lagcc_temp,'k','linewidth',1.5);
                %   plot(lags,lagcc_temp,'k','linewidth',1.5);
                title(['peak lag= ' num2str(lag) ]);
                grid on
                %   xlabel('listener precedes      tr (1.5s)     speaker precedes');
                
                subplot(3,2,3:4);
                plot(speaker_roi,'k','linewidth',1.5);
                title(['speaker ' rname]);
                grid on;
                xlabel('TR');
                
                subplot(3,2,5:6);
                
                plot(listener_froi,'k','linewidth',1.5);
                title({['listener ' frname ],['fROI centroid ['  num2str(centroid_mni(1)) ',' num2str(centroid_mni(2)) ',' num2str(centroid_mni(3)) ']' ]});
                grid on
                xlabel('TR');
                
                print(gcf,[expdir '/graph/' exp '/isfc_seed/' relation '/lagcorr_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_'  num2str(lag)  '_' rname '_' frname  '_['  num2str(centroid_mni(1)) ',' num2str(centroid_mni(2)) ',' num2str(centroid_mni(3)) ']' '_test.png' ],'-dpng');
                close gcf
            end
        end
    end
    
end




