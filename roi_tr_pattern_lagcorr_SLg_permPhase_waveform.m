
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;

rnames_selected={'vmPFC'};
for ei=[2];%2:4;
    exp=experiments{ei};
    
    for lagi=2%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_betaPeaks' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
        table(rnames,lag_neg,lag_pos,lag_both);
        
      %   load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'p_fdr');
        
        for ri=10;%1:length(rnames);
            rname=rnames{ri};
            
          %  if p_fdr(ri)<.05
          %        if ismember(rname,rnames_selected);
                %  if sig_betaClass(ri)==1;
                
                r_temp=r(ri,:);
                r_perm_temp=squeeze(r_perm(ri,:,:))';
                
                figure
                plot(lags,r_temp,'r');
                hold on;
                quantilePlot(r_perm_temp,lags,'k',0.3,[0.05/(length(rnames)*length(lags)) 1-0.05/(length(rnames)*length(lags))]);
                title({exp,rname})
                %    line([bPeakLags(ri,1) bPeakLags(ri,1)],[get(gca,'ylim')],'color','k');
                
                grid on
                hold off
                pause
                close gcf
                if ei==1;
                    xlabel('Speaker precedes-----TR(1.5s)-----Listeners precede','fontsize',10); ylabel('Beta values');
                end
         %   end
            
        end
    end
end

