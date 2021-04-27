close all
loc='mypc';
set_parameters
cols=jet(7);
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];
cols=cols.*0.95;
froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','DMNb'};
types={'wd','sn','pr'};
type_names={'word','sentence','paragraph'};
win_used=-5:20;

for ei=1:2;
    exp=exp_parameters.experiments{ei};
    
    %% by linguistic types
    fsize=[33 15];
    figure('unit','centimeter','position',[0 0 fsize]);
    
    load([expdir exp '/sound/onsets.mat']);
    
    for tpi=1:length(types);
        type=types{tpi};
        eval(['dur=dur_tr_' type ';']);
        [~, ri]=sort(dur);
        dur=round(dur);
        load(['Y:\claire\speaker-listener\' exp '\fmri\erp\network\' froidir '\erp_' type '.mat']);
        
        % y, Y, or Resid
        erps=erps_onset;
        
        % exclude shared boundaries and short boundaries
        if tpi==1;
            epi_excluded=(ismember(offsets_tr_wd,[offsets_tr_sn; offsets_tr_pr]) ) ;
            epi_excluded(dur<=0)=1;
        elseif tpi==2;
            epi_excluded=(ismember(offsets_tr_sn,offsets_tr_pr) );
            epi_excluded(dur<=1)=1;
        elseif tpi==3;
            % epi_excluded=(dur>=30);
            epi_excluded=[];
        end
        erps(:,:,:,epi_excluded)=NaN;
        
        %   dur_types=unique(dur);
        erps_bined=nan(size(erps));
        dur_types=unique(dur);
        erps_bined=erps_bined(:,:,:,1:length(dur_types));
        for di=1:length(dur_types);
            erps_bined(:,:,:,di)=nanmean(erps(:,:,:,dur==dur_types(di)),4);
        end
        erps=erps_bined;
        erps=erps(:,:,:,~isnan(squeeze(erps(1,1,1,:))));
        
        erps=erps(:,ismember(win,win_used),:,:);
        % baseline
        erps=erps-nanmean(erps,2);
        
        % average across subjects, CI across constituents
        erps=squeeze(nanmean(erps,3));
        cols=parula(size(erps,3));
        for ni=1:6;
            subplot(3,6,(tpi-1)*6+ni);
         %    temp=squeeze(zscore(erps(ni,:,:),0,2))';
          temp=squeeze(erps(ni,:,:))';
            if size(temp,2)==1; temp=temp'; end
            
            imagesc(temp);
            axis off
            %  imagesc(squeeze(erps(ni,:,:,:))');
            %             for di=1:size(erps,3);
            %                 plot(win_used,squeeze(erps(ni,:,di)),'color',cols(di,:),'linewidth',1.5);
            %                 hold on
            %             end
            %             hold off
            %             xlim([min(win_used) max(win_used)]);
            
            title(networks{ni});
        end
        axis on
        set(gca,'xtick',find(ismember(win_used,[-5:5:20])),'xticklabels',win_used(ismember(win_used,[-5:5:20])),'ytick',[]);
        ylabel({[type_names{tpi} 's'], 'ordered by length','long     short'})
        xlabel('time to onset (sec)');
    end
    
end
