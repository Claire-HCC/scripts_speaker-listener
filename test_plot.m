
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};


figure
for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
        subplot(2,4,ei);
        imagesc(zscore(r,0,2))
        title({exp,'temporal-spatial lagcorr'})
        set(gca,'xtick',1:length(lags),'xticklabel',lags);
        xlabel(['speaker precedes<---------->listeners precede']);
        
        subplot(2,4,ei+4)
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
        imagesc(zscore(r,0,2))
        title({exp,'temporal lagcorr'})
                set(gca,'xtick',1:length(lags),'xticklabel',lags);
        xlabel(['speaker precedes<---------->listeners precede']);
    end
end


figure
for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr');
        r_temp=r;
        r_temp(p_fdr>.05 | isnan(p_fdr))=NaN;
        subplot(2,4,ei);
        ris=find(nansum(r_temp,2)~=0);
        imagesc(r_temp)
        set(gca,'ytick',ris,'yticklabel',rnames(ris));
        set(gca,'xtick',1:length(lags),'xticklabel',lags);
        xlabel(['speaker precedes<---------->listeners precede']);
        title({exp,'temporal-spatial lagcorr'})
        

      load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr');
          r_temp=r;
        r_temp(p_fdr>.05 | isnan(p_fdr))=NaN;
        subplot(2,4,ei+4);
        ris=find(nansum(r_temp,2)~=0);
        imagesc(r_temp)
        set(gca,'ytick',ris,'yticklabel',rnames(ris));
        set(gca,'xtick',1:length(lags),'xticklabel',lags);
        xlabel(['speaker precedes<---------->listeners precede']);
        title({exp,'temporal lagcorr'})
    end
end

fsize=[20 20];
cols=jet(21);
for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','onsets');
        [onsets_sorted, ris_sorted]=sort(onsets,'ascend');
        rn=max(find(~isnan(onsets_sorted)));
        
        figure('unit','centimeter','position',[0 0 fsize]);
       
        hold on
        for ri=1:rn;
            
            plot(lags,r(ris_sorted(ri),:),'linewidth',2,'color',cols(onsets_sorted(ri)+11,:));
        end
        hold off
       legends = cellfun(@(x,y) [x ': ' num2str(y)],rnames(ris_sorted(1:rn)),mat2cell(onsets(ris_sorted(1:rn)),ones(rn,1),1),'UniformOutput',0);
        legend(legends,'location','bestoutside');
        legend boxoff
        title(exp)
        grid on
        set(gca,'fontsize',12);
         line([0 0],get(gca,'ylim'),'color','k')

    end
end