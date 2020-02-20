
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;
binSize=40;

for ei=3;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','onsets','onsetsOrder');
        [~,tn,~]=size(r);
        slope=nan([length(rnames) tn]);
        slope_p=nan([length(rnames) 1]);
        
        for ri=1:length(rnames);
            % restricted to speaker precedes
            for ti=1:tn;
                
                % ts=find(~isnan(onsets(ri,:)) & onsets(ri,:)<=0);
                ts=intersect(ti:binSize:tn,find(~isnan(onsets(ri,:)) & onsets(ri,:)<=0));
                % only in rois with more than 10 onsets
                if length(ts)>=5;
                    [slope(ri,ti), ~]=corr(ts',onsets(ri,ts)','type','spearman');
                end
            end
        end
        
        slopez=atanh(slope);
        [h,p,ci,stats] = ttest(slopez');
        roin=sum(nansum(~isnan(slope),2)>0);
        p_fdr=nan(size(p));
        [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
        
        
        
        save([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],
        
    end
end

ris=find(p_fdr<.05);
for i=1:length(ris);
    ri=ris(i);
    rname=rnames{ri};
        ts=intersect(1:tn,find(~isnan(onsets(ri,:)) & onsets(ri,:)<=0));
    figure;

    imagesc(squeeze(r(ri,:,:)))
    hold ;
    scatter(onsets(ri,:)+11,1:tn,'k','filled');
    title(sprintf('%s, R=%.2f, p-fdr=%.3f',rname,nanmean(slope(ri,:),2),p_fdr(ri)));
    hold off
    pause
    close gcf;
end


