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

for ei=[1 2 4 9:12];
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(rnames);
            seed=rnames{sdi};
            f=([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '.mat']);
            if exist(f);
                load(f,'r','keptT');
                if sdi==1;
                    rz=nan(size(r,1),size(r,1),size(r,2),size(r,3));
                end
                rz(sdi,:,:,:)=atanh(r);
            end
        end
        rzm=nanmean(rz,4);
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
        p=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>rzm,4);
        pfwe=p*length(p(:))/length(lags);
        
        peaks=nan(size(rzm,1),size(rzm,2));
        peakLags=nan(size(rzm,1),size(rzm,2));
        peaks_pfwe=nan(size(rzm,1),size(rzm,2));
        peakLags_pfwe=nan(size(rzm,1),size(rzm,2));
        peaks_p=nan(size(rzm,1),size(rzm,2));
        peakLags_p=nan(size(rzm,1),size(rzm,2));
        for sdi=1:length(rnames);
            for tgi=1:length(rnames);
                temp=squeeze(rzm(sdi,tgi,:));
                [pks, locs]=findpeaks(temp);
                locs=locs(pks>0);
                pks=pks(pks>0);
                if ~isempty(pks);
                    [~,loci]=(min(abs(locs-find(lags==0))));
                    peakLags(sdi,tgi)=lags(locs(loci));
                    peaks(sdi,tgi)=temp(locs(loci));
                    
                    % locs(~ismember(locs,find(squeeze(pfwe(sdi,tgi,:)<0.05))))=[];
                    % if ~isempty(locs);
                    %    [~,loci]=(min(abs(locs-find(lags==0))));
                    if p(sdi,tgi,locs(loci))<0.05;
                        peakLags_p(sdi,tgi)=lags(locs(loci));
                        peaks_p(sdi,tgi)=temp(locs(loci));
                    end
                    
                    if pfwe(sdi,tgi,locs(loci))<0.05;
                        peakLags_pfwe(sdi,tgi)=lags(locs(loci));
                        peaks_pfwe(sdi,tgi)=temp(locs(loci));
                    end
                end
                
            end
            
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksNearest' ],...
            'rnames','rzm','rz','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','peaks_p','peakLags_p');
        
    end
end

% %% find the peak nearest to lag 0 instead of the absolute peak
% clear all
% % loc='cluster';
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% lags_tested={-10:10,-15:15. -20:20};
% seed='aANG_R';
% 
% fsize=[20 20];
% figure('unit','centimeter','position',[0 0 fsize]);
% fi=1;
% for ei=[1 2 4 9:12];
%     exp=experiments{ei};
%     
%     for lagi=2;%1:length(lags_tested);
%         lags=lags_tested{lagi};
%         
%         load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksNearest' ],...
%             'rnames','rzm','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','peaks_p','peakLags_p');
%         nonnanN=sum(~isnan(peakLags_pfwe));
%         
%         if ei==1;
%             ri=find(strcmp(seed,rnames));
%             peakLags_0=squeeze(peakLags_pfwe(ri,:));
%             peakLagsM=nanmean(peakLags_pfwe);
%             [~,ris_sherlock]=sort(peakLagsM);
%             ris=ris_sherlock(nonnanN(ris_sherlock)>0);
%             
%             nii=roiTable2wholeBrainNii_mor([roi_ids(ris) peakLagsM(ris)']);
%             save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/peakLagsM.nii']);
%         end
%         ris=ris_sherlock(nonnanN(ris_sherlock)>0);
%         
%         imAlpha=ones(size(peakLags_pfwe(ris,ris)));
%          imAlpha(isnan(peakLags_pfwe(ris,ris)))=0;
% 
%         subplot(3,3,fi)
%         im=imagesc(peakLags_pfwe(ris,ris),[-5 5]); colormap jet
%                 set(im,'AlphaData',imAlpha)
%         % set(gca,'ytick',find(ris==ri),'yticklabels',seed,'xtick',[]);
%         set(gca,'ytick',[],'xtick',[]);
%         title([upper(exp(1)) strrep(exp,'_',' ')]);
%         set(gca,'fontsize',14,'color','k')
%         
%         fi=fi+1;
%     end
%     xlabel('target ROI')
%     ylabel('seed ROI')
%     
% end
% 
