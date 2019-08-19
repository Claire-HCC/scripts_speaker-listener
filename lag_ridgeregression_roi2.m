clear all
close all
set_parameters

lags=[-10:10]; %scan

% for ei=1%:2;
%     exp= experiments{ei};
%     epi_s=load([expdir experiments{ei} '/fmri/mat/roi/' exp '_speaker_rois.mat' ],'data');
%     epi_l=load([expdir experiments{ei} '/fmri/mat/roi/' exp '_listener_rois.mat' ],'data');
%
%     subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/listener*mat']));
%
%     roin=size(epi_s.data,1);
%     B_sl=zeros(roin,length(lags),length(subjects));
%     F_sl=zeros(roin,length(subjects));
%     p_sl=zeros(roin,length(subjects));
%     B_ll=zeros(roin,length(lags),length(subjects));
%     F_ll=zeros(roin,length(subjects));
%     p_ll=zeros(roin,length(subjects));
%
%     for si=1:length(subjects);
%         subj=subjects{si};
%
%         speaker=zscore(epi_s.data')';
%         listener_self=zscore(epi_l.data(:,:,si)')';
%         listener_others=zscore(mean(epi_l.data(:,:,:),3)')';
%         disp('Running coupling filter and lag correlation analysis...');
%
%         [B,F,p] = lag_ridgeregression3(listener_self(:,:), speaker(:,:), [lags(1) lags(end)]);
%         B_sl(:,:,si)=B;
%
%         [B,F,p]= lag_ridgeregression3(listener_self(:,:), listener_others(:,:), [lags(1) lags(end)]);
%         B_ll(:,:,si)=B;
%
%     end
%
%     clear B F p
%     B=B_sl; F=F_sl; p=p_sl;
%     save([expdir experiments{ei} '/fmri/mat/roi/isc/isc_sl.mat'],'B','F','p');
%     B=B_ll; F=F_ll; p=p_ll;
%     save([expdir experiments{ei} '/fmri/mat/roi/isc/isc_ll.mat'],'B','F','p');
% end


rois_selected={'vmPFC','dPCC','dPCUN','pIFG_L','dPreCG_L','HG_L','aOFC_L','DLPFC_R'};%,'aANG_L'; ,'PMC_L', ,'SMA_L', ,'IPL_R'
rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
r_inds=(find(contains(rtable.region,rois_selected)));
% r_inds=r_inds(1:6);

%
for ei=1%:2;
    exp= experiments{ei};
    load([expdir experiments{ei} '/fmri/mat/roi/isc/isc_sl.mat'],'B');
    B_sl=B;
    load([expdir experiments{ei} '/fmri/mat/roi/isc/isc_ll.mat'],'B');
    B_ll=B;
    
    
    for ri=1:2;%:length(r_inds);
        r_ind=r_inds(ri);
        rname= rtable.region(r_ind);
        
        B_sl_roi=squeeze(B_sl(r_ind,:,:))';
        B_ll_roi=squeeze(B_ll(r_ind,:,:))';
        
        B_ll_roi_ci=(ci(B_ll_roi,0.95));
        B_ll_roi_m=mean(B_ll_roi);
        
        B_sl_roi_ci=(ci( B_sl_roi,0.95));
        B_sl_roi_m=mean( B_sl_roi);
        
        figure;
        boundedline(lags,B_ll_roi_m,B_ll_roi_ci,'b',lags, B_sl_roi_m,B_sl_roi_ci,'r','transparency', 0.3)
        
        title(rname);
        legend({'otherListeners2OneListener','Speaker2OneListener'})
        xlabel(['speaker precedes    (scan)    listener precedes' ],'color','k');
    end
end

% 
% % individual plot
% for ei=1%:2;
%     exp= experiments{ei};
%     load([expdir experiments{ei} '/fmri/mat/roi/isc/isc_sl.mat'],'B');
%     B_sl=B;
%     load([expdir experiments{ei} '/fmri/mat/roi/isc/isc_ll.mat'],'B');
%     B_ll=B;
%     
%     for ri=1:2%1:length(r_inds);
%         r_ind=r_inds(ri);
%         rname= rtable.region(r_ind);
%         B_sl_roi=squeeze(B_sl(r_ind,:,:))';
%         B_ll_roi=squeeze(B_ll(r_ind,:,:))';
%         
%         figure;
%         
%         for si=1:size(B_ll,3);
%             subplot(4,5,si);
%             
%             B_ll_roi_ci=(ci(B_ll_roi,0.95));
%             B_ll_roi_m=mean(B_ll_roi);
%             boundedline(lags,B_ll_roi_m,B_ll_roi_ci,'b','transparency', 0.3)
%             
%             hold on
%             X=B_sl_roi(si,:);
%             plot(lags,X,'r');
%             %     ylim([-0.1 0.5])
%             title(['subject' num2str(si)]);
%             hold off
%         end
%         %    xlabel(['speaker precedes       listener precedes'
%         %    ],'color','r');
%         xlabel('scan');
%         subplot(4,5,si+1);
%         t=text([0 0 ],[0.5 0.7],{'otherListeners2OneListener','Speaker2OneListener'});
%         t(1).Color='b';
%         t(2).Color='r';
%         axis off;
%         title(rname);
%         
%     end
% end
% 
% 
