close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=2;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','peakLags_pfdr','peakLags');
        
        fsize=[35 28];
        
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags_pfdr);
        ris(isnan(tmp))=[];
        ris=ris(peakLags_pfdr(ris)>0);
        
        for i=1:length(ris);
            ri=ris(i);
            rname=rnames{ri};
            r_temp=r(ri,:);
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
        xlabel('Speaker precedes-----TR(1.5s)-----Listeners precede','fontsize',10);
        ylabel('R-value');
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',12)
    end
end

