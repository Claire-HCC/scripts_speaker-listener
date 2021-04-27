close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table','roi_newOrder');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};
seeds={'HG_L','pANG_L','vPCUN'};

rgb_surface=repmat(rgb,1,1,21);
rgb_surface=permute(rgb_surface,[1 3 2]);

surface(rzm(roi_newOrder,:),rgb_surface(roi_newOrder,:,1:3),'alphadata',p(roi_newOrder,:),'facealpha','flat','meshstyle','row','edgecolor','flat','linewidth',2)
surface(zeros(size(rzm)),'facealpha',0,'edgealpha',0.3)
surface(rzm(roi_newOrder,:),rgb_surface(roi_newOrder,:,1:3),'alphadata',(1-(p(roi_newOrder,:))+0.5)/1.5,'facealpha','flat','edgecolor','flat','meshstyle','none','alphadatamapping','none')


% fsize=[18 30];
% figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize],'color',[0.9 0.9 0.9]);
% sp=reshape(1:6,2,3)';

for ei=3:4;%:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=3;%1:length(seeds);
            seed=seeds{sdi};
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],'p','rnames','r','lags','peakLags_pfdr','peaks_pfdr','peakLags','peaks','peakLags_pfdr','peaks_pfdr','pfdr');
            
            rzm=nanmean(atanh(r),3);
            % to reverse the direction of information flow in figures.
            % rzm=fliplr(rzm);
            
            fsize=[9 9];
            figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize],'color',[0.9 0.9 0.9]);
            
            %     subplot(3,2,sp(sdi,ei-2));
            hold on;
            [tmp,ris]=sort(peakLags_pfdr);
            ris(isnan(tmp) |peaks_pfdr(ris)<0)=[];
            
            rgb=roi_table.rgb;
            %         rgb_temp=hot(length(lags));
            %         rgb=rgb_temp(nansum([peakLags repmat((length(lags)-1)/2,length(rnames),1)],2)+1,:);
            %         rgb(isnan(peakLags),:)=NaN;
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
       
          
             set(gca,'color','k','XColor','w','YColor','w','GridColor','w')  ;
             set(gca,'fontsize',16,'xtick',[min(lags):5:max(lags)]);
      %   title([upper(exp(1)) exp(2:end)],'color','w');
                title(strrep(seed,'_',' '),'color','w');
            grid on
            ylabel('R(z)','color','w');
          
            xticklabels([min(lags):5:max(lags)]*tr(ei))
            xlabel({['before seed--Lag(sec)--after seed']},'fontsize',14,'color','w');
        end
    end
end

