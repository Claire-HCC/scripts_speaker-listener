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
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','peakLags_pfdr','onsetsOrder','peak','peakLags_pfdr');
        r=mean(r,3);
        
        fsize=[20 18];
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
        [tmp,ris]=sort(peakLags_pfdr);
        ris(isnan(tmp))=[];
        for i=1:length(ris);
            ri=ris(i);
            rname=rnames{ri};
            r_temp=r(ri,:);
            plot(lags,r_temp,'linewidth',2,'color',cols(peakLags_pfdr(ri)+(length(lags)-1)/2+1,:));
            hold on
        end
        
        
        legends=cellfun(@(x,y) [x ' ' num2str(y)],rnames(ris),mat2cell(peakLags_pfdr(ris),ones(length(ris),1)),'Uniformoutput',0);
        legend(legends,'Location','eastoutside')
        legend boxoff
        
        colormap(cols);
        colorbar('location','southoutside','ticks',[]);
        
        set(gca,'fontsize',12)
    end
end

