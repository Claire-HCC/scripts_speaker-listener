close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10,  -40:40};
roi_ids=cell2mat(roi_table.id);
seeds={'HG_L','pANG_L','vPCUN'};

for ei=4%:4;
    exp=experiments{ei};
    
    for lagi=1;%2:3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        cols= jet(length(lags)+ceil(length(lags)/4)+1);
        cols = cols([1:length(lags)] + ceil(length(lags)/4)+1,:);
        
        for sdi=2:3;%1:length(seeds);
            seed=seeds{sdi};
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],'rnames','r','lags','peakLags_pfdr','peaks_pfdr','peakLags','peaks','peakLags_pfdr','peaks_pfdr','pfdr');
            
            rzm=nanmean(atanh(r),3);
            % to reverse the direction of information flow in figures.
            % rzm=fliplr(rzm);
            
            
            fsize=[9 9];
            figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize],'color',[0.9 0.9 0.9]);
            
            hold on;
            
            [tmp,ris]=sort(peakLags_pfdr);
            ris(isnan(tmp) |peaks_pfdr(ris)<0)=[];
            
            rgb=roi_table.rgb;
            %         rgb_temp=hot(length(lags));
            %         rgb=rgb_temp(nansum([peakLags repmat((length(lags)-1)/2,length(rnames),1)],2)+1,:);
            rgb(isnan(peakLags),:)=NaN;
            hold on;
            
            for ri=1:length(rnames);
                if ~ismember(ri,ris);
                    rname=rnames{ri};
                    r_temp=rzm(ri,:);
                    plot(lags,r_temp,'linewidth',1.5,'color',[0.8 0.8 0.8]);
                end
            end
            for i=1:length(ris);
                ri=ris(i);
                rname=rnames{ri};
                r_temp=rzm(ri,:);
                plot(lags,r_temp,'linewidth',2,'color',rgb(ri,:));
            end
            
            grid on
            
            set(gca,'color','k','XColor','w','YColor','w','GridColor','w')  ;
            set(gca,'fontsize',16,'xtick',[min(lags):5:max(lags)]);
            title(strrep(seed,'_',' '),'color','w');
            xticklabels([min(lags):5:max(lags)]*tr(ei))
            ylabel('R(z)','color','w');
            xlabel({['before seed--Lag(sec)--after seed']},'fontsize',14,'color','w');
        end
    end
end

