clear all
close all
set_parameters

lags=[-10:10]; %scan

rois_selected={'vmPFC','dPCC','dPCUN','pIFG_L','dPreCG_L','HG_L','aOFC_L','DLPFC_R'};%,'aANG_L'; ,'PMC_L', ,'SMA_L', ,'IPL_R'
rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
r_inds=(find(contains(rtable.region,rois_selected)));
% r_inds=r_inds(1:6);

figure;
for ei=1%:2;
    exp= experiments{ei};
    epi_s=load([expdir experiments{ei} '/fmri/mat/roi/' exp '_speaker_rois.mat' ],'data');
    epi_l=load([expdir experiments{ei} '/fmri/mat/roi/' exp '_listener_rois.mat' ],'data');
    
    segn=1;
    seg_len=size(epi_s.data(r_inds,:),2)/segn;
    for seg=1:segn;
        segt=((seg-1)*seg_len+1):(seg*seg_len);
        speaker=epi_s.data(r_inds,segt);
        listeners=epi_l.data(r_inds,segt,:);
        subjn=size(listeners,3);
        
        disp('Running coupling filter and lag correlation analysis...');
        
        [ B_regAll, B_regAuto, rss_regAll, rss_regAuto] = lag_ridgeregression2(speaker, listeners, [lags(1) lags(end)],2);
        gplot( B_regAll, B_regAuto, rss_regAll, rss_regAuto,rtable,r_inds,lags);
        
    end
end



%        individual plot

subplot_map=reshape(1:(4*length(r_inds)),length(r_inds),4);
subplot_map=subplot_map';

for subj=1:size(B_regAll,3);
    if mod(subj,4)==1; figure; end
    
    for ri=1:length(r_inds);
        
        r_ind=r_inds(ri);
        
        si=mod(subj,4); if si==0; si=4; end
        subplot(4,length(r_inds),subplot_map(si,ri));
        %    subplot(segn,length(r_inds),ri+length(r_inds)*(seg-1));
        
        plot(lags(lags~=0),B_regAll(ri,:,subj));
        hold on
        plot(lags(lags~=0),B_regAuto(ri,:,subj));
        
        %  legend({['RSS=' num2str(rss_regAll(ri,subj))],['RSS=' num2str(rss_regAuto(ri,:,subj))]},'Location','south');
        
        grid on
        hold off
        
        title(rtable.region{r_ind});
        
        if ri==1;
            %  legend({['speaker2listener, RSS=' num2str(round(mean(rss_regAll(ri,:)))) ],['listener2listener, RSS=' num2str(round(mean(min(rss_regAuto(ri,:,:)))))]},'Location','south');
            ylabel({['Subj' num2str(subj)],'coef'});
            
        end
        
    end
end
legend({['listener2speaker' ],['speaker2speaker']},'Location','south');


% group plot
function gplot( B_regAll, B_regAuto, rss_regAll, rss_regAuto,rtable,r_inds,lags);
figure;


for ri=1:size( B_regAll);
    r_ind=r_inds(ri);
    
    subplot(2,length(r_inds)/2,ri);
    
    B_regAll_m=mean(B_regAll(ri,:,:),3);
    B_regAll_ci=(ci(squeeze(B_regAll(ri,:,:))',0.95));
    confplot(lags(lags~=0),B_regAll_m,B_regAll_ci);
    
    hold on
    B_regAuto_m=mean(B_regAuto(ri,:,:),3);
    plot(lags(lags~=0),B_regAuto_m)
    
    
  %  legend({['RSS=' num2str(round(mean(rss_regAll(ri,:))))],['RSS=' num2str(round(mean(min(rss_regAuto(ri,:,:)))))]},'Location','south');
    
    title(rtable.region{r_ind});
    grid on
    hold off
    
    if ri==1;
        %    legend({['multiple regressor, RSS=' num2str(round(mean(rss_regAll(ri,:)))) ],['one regressor, RSS=' num2str(round(mean(min(rss_regAuto(ri,:,:)))))]},'Location','south');
        %            ylabel({['story segment ' num2str(seg)],'regression coef'});
        %   xlabel('lags (scan)');
    end
    
end
legend({['listener2speaker' ],['speaker2speaker']},'Location','south');
end


