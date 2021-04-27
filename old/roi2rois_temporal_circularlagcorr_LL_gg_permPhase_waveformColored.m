close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

seed='vPCUN';
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','peakLags_pfdr','peakLags','peaks_pfdr','peakLags_pfwe','peaks_pfwe');
        rzg=nanmean(atanh(r),3);
        fsize=[30 25];
        
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags_pfwe);
        ris(isnan(tmp) | peaks_pfwe(ris)<0)=[];
        
        for i=1:length(ris);
            ri=ris(i);
           
                rname=rnames{ri};
                r_temp=rzg(ri,:);
                subplot(1,3,[1 2 ])
                plot(lags,r_temp,'linewidth',2,'color',cols(peakLags_pfwe(ri)+(length(lags)-1)/2+1,:));
                hold on
                
                subplot(1,3,3);
                if i==1 | sign(peakLags_pfwe(ri))~= sign(peakLags_pfwe(ris(i-1))); text_y=30 ;
                else text_y=text_y-1 ;
                end
                xlim([-1 1.5])
                ylim([0 30])
                text(sign(peakLags_pfwe(ri)),text_y,'-','fontsize',16,'fontweight','bold','HorizontalAlignment','left','VerticalAlignment','middle','color',cols(peakLags_pfwe(ri)+(length(lags)-1)/2+1,:));
                text(sign(peakLags_pfwe(ri))+0.1,text_y,[ strrep(rname,'_',' ') ': ' num2str(peakLags_pfwe(ri))  ],'color','k','HorizontalAlignment','left','fontsize',12);

        end
        axis off
        
        subplot(1,3,[1 2 ])
        title({exp})
        grid on
        xlabel({['later than ' seed '---TR---earlier than ' seed],['<-----Information Flow <--------']},'fontsize',10);
        ylabel('R(z)');
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',14)
    end
end

