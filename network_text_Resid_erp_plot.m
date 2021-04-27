clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6_audPausesResid';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};


cols=jet(7);
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];
cols=cols.*0.95;

types={'wd','sn','pr'};
type_names={'word','sentence','paragraph'};
eis=[1 2 6];


for eii=1;%length(eis);
    ei=eis(eii);
    exp=exp_parameters.experiments{ei};
    
    fsize=[27 27]
    figure('unit','centimeter','position',[0 0 fsize]);
    
    spi=1;
    for tpi=1:3;%:length(types);
        type=types{tpi};

        load(['Y:\claire\speaker-listener\' exp '\fmri\erp\network\' froidir '\erp_' type '.mat']);;
        erps_cell={erpsy_offset,erpsY_offset,erps_offset};
        win_used=-30:30;
        
        for epi=1:length(erps_cell);
            ax(spi)=subplot(length(types),length(erps_cell),spi)
            erps=erps_cell{epi};
            erps=erps(:,ismember(win,win_used),:,:);
            erps=erps-nanmean(erps,2);
            
            % average across subject, CI across boudnaries
            erps=squeeze(nanmean(erps,3));
            
            % ciplot_claire(squeeze(erps_aud(1,:,:))',win,[0.3 0.3 0.3] ,0.1);
            % hold on
            for sdi=1:length(networks);
           ciplot_claire(squeeze(erps(sdi,:,:))',win_used,[cols(sdi,:)] ,0.1);
                hold on;
            end
            set(gca,'fontsize',12)
            grid on;
            xlim([min(win) max(win)]);
            xlim([min(win_used) max(win_used)]);
            %             if epi==1;
            %                 yl=get(gca,'ylim');
            %             end
      
            set(gca,'xtick',-40:10:40,'xticklabels',[-40:10:40]*1.5);
            xlabel({'Time relstive to',[type_names{tpi} ' boundary (sec)']});
            ylabel('fMRI signal');
            title([upper(exp(1)) strrep(exp(2:end),'_',' ')]);
            
            spi=spi+1;
        end
        
    end
         linkaxes(ax);
end
