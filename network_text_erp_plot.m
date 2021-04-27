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
win_used=-5:15;

for ei=[6];
    exp=exp_parameters.experiments{ei};
    
    %% by linguistic types
    fsize=[28 9];
    figure('unit','centimeter','position',[0 0 fsize]);
    
    load([expdir exp '/sound/onsets.mat']);
    
    for tpi=3;%1:length(types);
        type=types{tpi};
        
        load(['Y:\claire\speaker-listener\' exp '\fmri\erp\network\' froidir '\erp_' type '.mat']);
        
        % y, Y, or Resid
        erps=erps_offset;
        
        % exclude shared boundaries and short boundaries
        if tpi==1;
            epi=(ismember(offsets_tr_wd,offsets_tr_sn) | ismember(offsets_tr_wd,offsets_tr_pr) ) ;
        elseif tpi==2;
            epi=(ismember(offsets_tr_sn,offsets_tr_pr) );
        elseif tpi==3;
            epi=[];
        end
        erps(:,:,:,epi)=NaN;
        erps=erps(:,:,:,~isnan(squeeze(erps(1,1,1,:))));
        
        erps=erps(:,ismember(win,win_used),:,:);
        
        % baseline
        % erps=erps-nanmean(erps(:,1:5,:,:),2);
        erps=erps-nanmean(erps,2);
        
        % average across subjects, CI across constituents
        %erps=squeeze(nanmean(erps,3));
        
        subplot(1,3,tpi);
        hold off;
        % ciplot_claire(squeeze(erps_aud(1,:,:))',win,[0.3 0.3 0.3] ,0.1);
        % hold on
        for sdi=1:length(networks);
            ciplot_claire(squeeze(erps(sdi,:,:))',win_used,[cols(sdi,:)] ,0.1);
            hold on;
        end
        
        xlim([min(win_used) max(win_used)]);
        set(gca,'xtick',-16:4:16,'xticklabels',[-16:4:16]*1.5);
        set(gca,'fontsize',14);
        grid on;
        xlabel('time to boundary (sec)');
        ylabel('fMRI signal');
        title(type_names{tpi});
    end
    ylim([-30 30])
    yl=get(gca,'ylim');
    for tpi=1:length(types);
        subplot(1,3,tpi);
        ylim(yl);
    end
    
end
