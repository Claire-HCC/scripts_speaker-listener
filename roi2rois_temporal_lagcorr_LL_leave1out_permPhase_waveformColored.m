close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

seed='vPCUN'
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','rz','lags','peakLags_pfdr','peakLags');
        rz=nanmean(rz,3);
        fsize=[35 28];
        
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags_pfdr);
        ris(isnan(tmp))=[];
        
        for i=1:length(ris);
            ri=ris(i);
            rname=rnames{ri};
            r_temp=rz(ri,:);
            subplot(1,4,[1 2 3])
            plot(lags,r_temp,'linewidth',2,'color',cols(peakLags_pfdr(ri)+(length(lags)-1)/2+1,:));
            hold on
            
            subplot(1,4,4);
            if i==1 | sign(peakLags_pfdr(ri))~= sign(peakLags_pfdr(ris(i-1))); text_y=40 ;
            else text_y=text_y-1 ;
            end
            xlim([-1 1.5])
            ylim([0 40])
            text(sign(peakLags_pfdr(ri)),text_y,'-','fontsize',20,'HorizontalAlignment','left','VerticalAlignment','middle','color',cols(peakLags_pfdr(ri)+(length(lags)-1)/2+1,:));
            text(sign(peakLags_pfdr(ri))+0.1,text_y,[ rname ': ' num2str(peakLags_pfdr(ri))  ],'color','k','HorizontalAlignment','left');
            
        end
        axis off
        
        subplot(1,4,[1 2 3])
        title({exp})
        grid on
        xlabel(['later than ' seed '---TR---earlier than ' seed],'fontsize',10);
        ylabel('R(z)');
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',14)
    end
end

