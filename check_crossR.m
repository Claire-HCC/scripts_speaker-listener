clear all
close all
set_parameters
rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
rois_selected={'HG_L','HG_R'};
rids=rtable.id(find(contains(rtable.region,rois_selected)));
lag=50;

for ei=1:2;
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/'  exp '_speaker_subj01.mat']);
    epi_s=zscore(data')';
    
    subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' exp '_listener*mat']));
    
    subploti=reshape(1:(length(subjects)*length(rids)),length(rids),length(subjects))';
    fsize=[10 30];
    h=figure('position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize,'unit','centimeter');
    for si=1:length(subjects);
        subj=subjects{si};
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi_l=zscore(data')';
        
        for ri=1:length(rids);
            rid=rids{ri};
            %             load([expdir '/roi_mask/mat/HG_L.mat']);
            %             rmask=data;
            load([expdir '/roi_mask/mat/' rid '.mat']);
            rmask=data;
            epi_s_masked=mean(epi_s(rmask,:))';
            epi_l_masked=mean(epi_l(rmask,:))';
            
            epi_s_masked=despike(epi_s_masked,[-3 3],5);
             epi_l_masked=despike(epi_l_masked,[-3 3],5);
            
            voln=min(length( epi_s_masked),length( epi_l_masked));
            %    figure; subplot(2,1,1); plot(zscore(aud)); xlim([0 voln]); subplot(2,1,2); plot(zscore(epi_masked));  xlim([0 voln]);
            lag=50;
            [r,lags]=xcorr(epi_s_masked(1:voln),epi_l_masked(1:voln),lag,'coeff');
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
    print(gcf,[expdir '/graph/fMRI_sanityCheck/' experiments{ei}  '_crossR.png'],'-dpng');
end


