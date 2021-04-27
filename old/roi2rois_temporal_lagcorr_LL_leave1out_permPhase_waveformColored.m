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
for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],'pfdr','rnames','r','lags','peakLags_pfdr','peakLags','peaks_pfdr','peakLags_pfdr','peaks_pfdr');
        rzg=nanmean(atanh(r),3);
        fsize=[30 24];
        
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags_pfdr);
        ris(isnan(tmp) |peaks_pfdr(ris)<0)=[];
        
        for i=1:length(ris);
            ri=ris(i);
            
            rname=rnames{ri};
            r_temp=rzg(ri,:);
            subplot(1,3,[1 2 ])
            plot(lags,r_temp,'linewidth',2,'color',cols(peakLags_pfdr(ri)+(length(lags)-1)/2+1,:));
            hold on
            
            subplot(1,3,3);
            if i==1 | sign(peakLags_pfdr(ri))~= sign(peakLags_pfdr(ris(i-1))); text_y=30 ;
            else text_y=text_y-1 ;
            end
            xlim([-0.5 1.5])
            ylim([0 30])
            text(sign(peakLags_pfdr(ri)),text_y,'-','fontsize',20,'fontweight','bold','HorizontalAlignment','left','VerticalAlignment','middle','color',cols(peakLags_pfdr(ri)+(length(lags)-1)/2+1,:));
            text(sign(peakLags_pfdr(ri))+0.1,text_y,[ strrep(rname,'_',' ') ': ' num2str(peakLags_pfdr(ri))  ],'color','k','HorizontalAlignment','left','fontsize',14);
            
        end
        axis off
        
        subplot(1,3,[1 2 ])
        title([upper(exp(1)) exp(2:end)])
        grid on
        
        ylabel('R(z)');
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',20)
        
        xlabel({['later than ' seed '---TR---earlier than ' seed],['<-----Information Flow <--------']},'fontsize',18);
    end
end

