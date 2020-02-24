close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_permSL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','peakLags2_pfdr','peakLags2');
        
        fsize=[20 25];
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags2);
        ris(isnan(tmp))=[];
        for i=1:length(ris);
            ri=ris(i);
            rname=rnames{ri};
            r_temp=r(ri,:);
            plot(lags,r_temp,'linewidth',2,'color',cols(peakLags2(ri)+(length(lags)-1)/2+1,:));
            hold on
        end
        title({exp})
        grid on
        xlabel('Speaker precedes-----TR(1.5s)-----Listeners precede','fontsize',10);
        ylabel('R-value');
        
        legends=cellfun(@(x,y) [x ' ' num2str(y)],rnames(ris),mat2cell(peakLags2(ris),ones(length(ris),1)),'Uniformoutput',0);
        legend(legends,'Location','eastoutside')
        legend boxoff
        
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',12)
    end
end

