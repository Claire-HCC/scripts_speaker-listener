
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=10;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname='dPCC';
ri=find(ismember(rnames,rname));

for ei=3%1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
    r2_ll=nanmean(r2,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL.mat'],'r2');
    r2_sl_perm=r2;
    r2_sl_perm_m=nanmean(r2,3);
    
    keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
     keptT2=keptT+binSize/2;
     
    r2_ll_temp=r2_ll(ri,keptT);
    r2_sl_temp=r2_sl(ri,keptT);
    r2_sl_perm_m_temp=r2_sl_perm_m(ri,keptT);
    
    cols=zeros(length(keptT),3);
    colx=r2_sl_temp> r2_sl_perm_m_temp;%((max(r2_sl_temp)-min(r2_sl_temp))/2+min(r2_sl_temp));
    coly=r2_ll_temp>median(r2_ll_temp);%((max(r2_ll_temp)-min(r2_ll_temp))/2+min(r2_ll_temp));
    
    for ti=1:length(keptT2);
        if  colx(ti)==1 & coly(ti)==1;
            cols(ti,:)=[1 0.7 0];
        elseif colx(ti)==1 & coly(ti)==0;
            cols(ti,:)=[0 1 0];
        elseif colx(ti)==0 & coly(ti)==0;
            cols(ti,:)=[0 0 1];
        elseif colx(ti)==0 & coly(ti)==1;
            cols(ti,:)=[0.7 0.1 1];
        end
    end
    
    % so that the line is in the center of the bin!!!
   
    fsize=[30 9];
    f=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    hold on
    for ti=1:length(keptT2);
        t=keptT2(ti);
        line([t t],[0 1],'color',[cols(ti,:) 0.3])
    end
    
    plot(keptT2,[r2_sl_temp],'r','linewidth',2);
    plot(keptT2,[r2_ll_temp],'b','linewidth',2);
    
    hold off
    title({[upper(exp(1)) exp(2:end) ', ' rnames{ri}]});%,['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
    xlim([min(keptT2) max(keptT2)])
    ylim([0 max(r2_sl_temp(:))+0.05]);
    
    xlabel('Time (TR)');
    %  ylabel('Neural Coupling');
    set(gca,'fontsize',14);
    h = findobj(gca,'Type','line');
    legend(h,'Listener-Listener Coupling','Speaker-Listener Coupling','location','northwest');
    legend boxoff
    set(gca,'ytick',[])
    set(gca,'ycolor','w');
    
    [aud, fs_aud] = audioread([expdir experiments{ei} '/sound/' exp '_listener_cropped.wav']);
    aud([min(keptT)*tr(ei)*fs_aud  (max(keptT)+binSize-1)*tr(ei)*fs_aud]);
    fs_vid =2;% frame rate must be larger than 1, so i divid 1.5 sec (one tr) in to 3 frames
    aud_perFrame = fs_aud / fs_vid;
 %  videoFWriter = vision.VideoFileWriter([expdir '/' exp '/fmri/animation/' rname '_pattern_binSize' num2str(binSize) '_herdCategories.avi'], 'FrameRate',fs_vid , 'AudioInputPort', true);
     videoFWriter = vision.VideoFileWriter([expdir '/' rname '_pattern_binSize' num2str(binSize) '_herdCategories.avi'], 'FrameRate',fs_vid , 'AudioInputPort', true,'FileFormat','AVI');

    for ti=1:20;%length(keptT2);
        t=keptT2(ti);
        if ti>1;
            delete(l);
        end
        l=line([t t],[0 1],'color','k');
        
        F=getframe(f);
        F=F.cdata;
        
        for tt=1:3;
            audind_s= floor(((keptT(ti) - 1)*(tr(ei)/0.5)+(tt-1) )*aud_perFrame ) + 1;
            audind_e= ceil(audind_s + aud_perFrame - 1);
            step(videoFWriter,F, aud(audind_s:audind_e));
        end
    end
    release(videoFWriter);
end


