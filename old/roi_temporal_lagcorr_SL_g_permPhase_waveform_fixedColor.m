close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','peakLags_pfdr','peaks_pfdr','peakLags','peaks','peakLags_pfdr','peaks_pfdr','pfdr');
        
        fsize=[15 16];
        figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize],'color',[0.9 0.9 0.9]);
        
        [tmp,ris]=sort(peakLags_pfdr);
        ris(isnan(tmp) |peaks_pfdr(ris)<0)=[];
        %         [tmp,ris]=sort(peakLags);
        %         ris(isnan(tmp) |peaks_pfdr(ris)<0)=[];
        %
        %roi_table.rgb;
        rgb_temp=hot(length(lags));
        rgb=rgb_temp(nansum([peakLags repmat((length(lags)-1)/2,length(rnames),1)],2)+1,:);
        rgb(isnan(peakLags),:)=NaN;
        hold on;
        
        
        % rgb=roi_table.rgb;
        rgb_temp=hot(length(lags));
        rgb=rgb_temp(nansum([peakLags repmat((length(lags)-1)/2,length(rnames),1)],2)+1,:);
        rgb(isnan(peakLags),:)=NaN;
        hold on;
        
        for ri=1:length(rnames);
            if ~ismember(ri,ris);
                rname=rnames{ri};
                r_temp=nanmean(atanh(r(ri,:,:)),3);
                plot(lags,r_temp,'linewidth',2,'color',[0.8 0.8 0.8]);
            end
        end
        for i=1:length(ris);
            ri=ris(i);
            rname=rnames{ri};
            r_temp=nanmean(atanh(r(ri,:,:)),3);
            plot(lags,r_temp,'linewidth',3,'color',rgb(ri,:));
        end
        title([upper(exp(1)) exp(2:end)])
        grid on
        ylabel('R(z)');
        set(gca,'fontsize',18)
        xlabel({['before speaker---Peak Lag (TR)---after speaker'],[ '<---information flow<---']},'fontsize',16);
    end
end

