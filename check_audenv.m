clear all
close all
set_parameters
rtable=readtable([expdir 'roi_mask/roi_id_region.txt']);
rois_selected={'HG_L','HG_R'};
rids=rtable.id(find(contains(rtable.region,rois_selected)));
lag=50;

for ei=1:2;
    load([expdir experiments{ei} '/sound/' experiments{ei} '_audenv.mat' ]);
    subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*mat']));
    
    subploti=reshape(1:(length(subjects)*length(rids)),length(rids),length(subjects))';
    fsize=[10 30];
    h=figure('position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize,'unit','centimeter');
    for si=1:length(subjects);
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi=data;
        
        for ri=1:length(rids);
            rid=rids(ri);
            load([expdir '/roi_mask/mat/roi' num2str(rid) '.mat']);
            rmask=data;
            epi_masked=mean(epi(rmask,:))';
            
            voln=min(length(epi_masked),length(aud));
            %    figure; subplot(2,1,1); plot(zscore(aud)); xlim([0 voln]); subplot(2,1,2); plot(zscore(epi_masked));  xlim([0 voln]);
            lags=-50:50;
            r=lagcorr(zscore(epi_masked(1:voln))',zscore(aud(1:voln))',lags);
            pt=lags(r==max(r))
            set(0, 'currentfigure', h);
            subplot(length(subjects),length(rids), subploti(si,ri));
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
    print(gcf,[expdir '/graph/fMRI_sanityCheck/' experiments{ei}  '_audenv.png'],'-dpng');
end


