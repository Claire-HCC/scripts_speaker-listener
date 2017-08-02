clear all
close all
set_parameters

lag=4;
tic

for ei=1%:2;
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/'  exp '_speaker_subj1.mat']);
    epi_s=zscore(data')';
    
    subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' exp '_listener*mat']));
    %  fsize=[10 30];
    %  h=figure('position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize,'unit','centimeter');
    for si=1:length(subjects);
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi_l=zscore(data')';
        
        
        mask=find(sum(epi_s')~=0 & sum(epi_l')~=0);
        
        for v = mask;
            [r,lags]=xcorr(epi_s(v,:),epi_l(v,:),lag,'coeff');
            pk=max(r);
            pt=lags(r==pk);
            isc_sl_peakR(v,1,si)=pk;
            isc_sl_peakT(v,1,si)=pt;
        end
    end
end

save([expdir exp '/fmri/mat/wholeBrain/isc_sl_peakR.mat' ],'isc_sl_peakR');
save([expdir exp '/fmri/mat/wholeBrain/isc_sl_peakT.mat' ],'isc_sl_peakT');
toc
