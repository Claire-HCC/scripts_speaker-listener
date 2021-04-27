%% ttest across subject
clear all
% loc='cluster';
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

for ir=1:length(interROItypes);
    interROItype=interROItypes{ir};
    
    if strcmp(interROItype,'roi');
        seeds_={''};
    else
        seeds_= cellfun(@(x) [x '_'],seeds,'UniformOutput',0);
    end
    
    for is=1:length(interSubjTypes);
        interSubjType=interSubjTypes{is};
        
        for ei=3:4;
            exp=experiments{ei};
            
            for sdi=1:length(seeds_);
                seed_=[seeds_{sdi}];
                
                for lagi=1%:length(lags_tested);
                    lags=lags_tested{lagi};
                    
                    myVideo = VideoWriter([expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr.avi']);
                    myVideo.FrameRate =tr(ei);
                    myVideo.Quality =100;
                    open(myVideo);
                    
                    for li=1:length(lags);
                        lag=lags(li);
                        
                        f=figure;
                        fs{1}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_lh_l.tif'];
                        fs{2}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_lh_m.tif'];
                        fs{3}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_rh_l.tif'];
                        fs{4}=[expdir exp '\fmri\temporal\lagcorr\tr\' interROItype '\' froidir '\' interSubjType '\' seed_ 'lag' num2str(min(lags)) '-' num2str(max(lags))  '_R_pfdr_lag' num2str(lag) '_rh_m.tif'];
                        montage(fs,'size',[4 1]);
                        set(gcf,'color','k');
                        text(0.17*mean(get(gca,'xlim')),mean(get(gca,'ylim')),'Lag','horizontalalignment','left','color','w','fontsize',20);
                        text(1.75*mean(get(gca,'xlim')),mean(get(gca,'ylim')),sprintf('%3.1f sec',lag*tr(ei)),'color','w','horizontalalignment','right','fontsize',20);
                        
                        F(li) = getframe(f);
                    
                        close gcf
                        
                    end
                    
                    writeVideo(myVideo,F);
                    close(myVideo)
                end
            end
        end
    end
end


