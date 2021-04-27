%% ttest across subject
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};
interROItypes={'roi','roi2rois'};
seeds={'vPCUN','pANG_L','HG_L'};
interSubjTypes={'SL_each','LL_leave1out'};

for ir=2;%1:length(interROItypes);
    interROItype=interROItypes{ir};
    
    if strcmp(interROItype,'roi');
        seeds_={''};
    else
        seeds_= cellfun(@(x) [x '_'],seeds,'UniformOutput',0);
    end
    
    for is=1;%1:length(interSubjTypes);
        interSubjType=interSubjTypes{is};
        
        for ei=3%:4;
            exp=experiments{ei};
            
            for sdi=1;%:length(seeds_);
                seed_=[seeds_{sdi}];
                
                for lagi=1%:length(lags_tested);
                    lags=lags_tested{lagi};
                    
                  
                   lagsShown=[-4 -2 0 2 4];
                   lagsShownN=length(lagsShown);
                    for i=1:lagsShownN;%1:length(lags);
                        lag=lagsShown(i);
                        
                        fs{i,1}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_lh_l.tif'];
                        fs{i,2}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_lh_m.tif'];
                       
                    end
                         figure
                    montage(fs(:),'size',[2 lagsShownN]);
                    x=min(get(gca,'xlim')):diff(get(gca,'xlim'))/(2*lagsShownN):max(get(gca,'xlim'));
                    x=x(2:2:end);
                    y=repmat(mean(get(gca,'ylim')),1,lagsShownN);
                 
                    text(x,y,cellfun(@(x) [x ' sec'],cellstr(num2str(lagsShown'*tr(ei))),'UniformOutput',0),'color','w','horizontalalignment','center','fontsize',20)

                    
                end
            end
        end
    end
end




