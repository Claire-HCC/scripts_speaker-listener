clear all
close all
set_parameters

%rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
rois_selected={'HG_L','HG_R'};

% rids=rtable.id(find(contains(rtable.region,rois_selected)));
lag=50;

type='env';
role='listener';

for ei=3:4%:2;
   
    subjects=cellstr(ls([expdir experiments{ei} '/fmri/timeseries/wholeBrain/' cropped '/*mat']));
    
    subploti=reshape(1:(length(subjects)*length(rids)),length(subjects),length(rids))';
    fsize=[30 10];
    h=figure('position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize,'unit','centimeter');
    
    for si=1:length(subjects);
        subj=subjects{si};
        
         load([expdir experiments{ei} '/sound/' experiments{ei} '_listener_cropped_aud' type '.mat' ]);
        
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' cropped '/' subj]);
        epi=data;
        
        for ri=1:length(rids);
            rid=rids{ri};
            load([expdir '/roi_mask/mat/' rid '.mat']);
            rmask=data;
            epi_masked=mean(epi(rmask,:))';
            
            voln=min(length(epi_masked),length(aud));
            %    figure; subplot(2,1,1); plot(zscore(aud)); xlim([0 voln]); subplot(2,1,2); plot(zscore(epi_masked));  xlim([0 voln]);
            lag=50;
           [r,lags]=xcorr(zscore(epi_masked(1:voln))',zscore(aud(1:voln))',lag,'coeff');
            pt=lags(r==max(r))
            set(0, 'currentfigure', h);
            subplot(length(rids),length(subjects), subploti(ri,si));
            plot(lags,r);
            grid on;
            line([pt pt],[-1 1]);
            ylim([-0.5 0.5]);
            
            if si==1;
                title({rois_selected{ri},strrep(subj,'_',' '),['peak at ' num2str(pt) ' scan']});
            else
                title({strrep(subj,'_',' '),['peak at ' num2str(pt) ' scan']});
            end
            hold on
        end
    end
    xlabel('scan');
 print(gcf,[expdir '/graph/fMRI_sanityCheck/' experiments{ei}  '_' cropped '_aud' type '.png'],'-dpng');
end


