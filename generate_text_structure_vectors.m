clear all
set_parameters;

for ei=1:2;
    exp=exp_parameters.experiments{ei};
    
    tb=readtable([expdir exp '/sound/transcripts_'  exp '_Claire.xlsx']);
    tb(strcmp(tb.word,'um'),:)=[];
    tb(strcmp(tb.word,'oh'),:)=[];
    tb(strcmp(tb.word,'ya'),:)=[];
    tb(strcmp(tb.word,'yeah'),:)=[];
    
    tb.tmax_msec=tb.tmax*1000;
    tb.tmin_msec=tb.tmin*1000;
    
    onsets_msec_wd=tb.tmin_msec;
    onsets_msec_sn=tb.tmin_msec(find([1; abs(diff(tb.sentencei))]>0  & tb.sentencei>0));
    onsets_msec_pr=tb.tmin_msec(find(([1; abs(diff(tb.paragraphi))]>0 |isnan([1; abs(diff(tb.paragraphi))]) & tb.paragraphi>0)));
    
    offsets_msec_wd=tb.tmax_msec;
    offsets_msec_sn=tb.tmax_msec(find([abs(diff(tb.sentencei)); 1]>0 & tb.sentencei>0));
    offsets_msec_pr=tb.tmax_msec(find(([abs(diff(tb.paragraphi)); 1]>0  | isnan([abs(diff(tb.paragraphi)); 1])& tb.paragraphi>0)));
    
    onsetsVectorMsec_wd=zeros(1,round(max(tb.tmax_msec)));
    for wi=1:length(onsets_msec_wd);
        
        onset=round(onsets_msec_wd(wi));
        offset=round(offsets_msec_wd(wi));
        
        onsetsVectorMsec_wd(onset)=1;
    end
    
    onsetsVectorMsec_sn=zeros(1,round(max(tb.tmax_msec)));
    for sni=1:length(onsets_msec_sn);
        onset=round(onsets_msec_sn(sni));
        offset=round(offsets_msec_sn(sni));
        onsetsVectorMsec_sn(onset)=1;
    end
    
    onsetsVectorMsec_pr=zeros(1,round(max(tb.tmax_msec)));
    for pri=1:length(onsets_msec_pr);
        
        onset=round(onsets_msec_pr(pri));
        offset=round(offsets_msec_pr(pri));
        onsetsVectorMsec_pr(onset)=1;
    end
    
    onsets_tr_wd=tb.tmin_tr;
    onsets_tr_sn=tb.tmin_tr(find([1; abs(diff(tb.sentencei))]>0 & tb.sentencei>0));
    onsets_tr_pr=tb.tmin_tr(find(([1; abs(diff(tb.paragraphi))]>0 |isnan([1; abs(diff(tb.paragraphi))]) & tb.paragraphi>0)));
    
    offsets_tr_wd=tb.tmax_tr;
    offsets_tr_sn=tb.tmax_tr(find([abs(diff(tb.sentencei)); 1]>0 & tb.sentencei>0));
    offsets_tr_pr=tb.tmax_tr(find(([abs(diff(tb.paragraphi)); 1]>0  | isnan([1; abs(diff(tb.paragraphi))])& tb.paragraphi>0)));
    
    onsetsVectorTr_wd=zeros(1,round(max(tb.tmax_tr)));
    for wi=1:length(onsets_tr_wd);
        
        onset=round(onsets_tr_wd(wi));
        offset=round(offsets_tr_wd(wi));
        
        if onset>0 & onset < length(onsetsVectorTr_wd);
            onsetsVectorTr_wd(onset)=1;
        end
    end
    
    onsetsVectorTr_sn=zeros(1,round(max(tb.tmax_tr)));
    for sni=1:length(onsets_tr_sn);
        onset=round(onsets_tr_sn(sni));
        offset=round(offsets_tr_sn(sni));
        if onset>0 & onset < length(onsetsVectorTr_wd);
            onsetsVectorTr_sn(onset)=1;
        end
    end
    
    onsetsVectorTr_pr=zeros(1,round(max(tb.tmax_tr)));
    for pri=1:length(onsets_tr_pr);
        
        onset=round(onsets_tr_pr(pri));
        offset=round(offsets_tr_pr(pri));
        if onset>0 & onset < length(onsetsVectorTr_wd);
            onsetsVectorTr_pr(onset)=1;
        end
    end
    dur_tr_wd=offsets_tr_wd-onsets_tr_wd;
    dur_tr_sn=offsets_tr_sn-onsets_tr_sn;
    dur_tr_pr=offsets_tr_pr-onsets_tr_pr;
    
    dur_msec_wd=offsets_msec_wd-onsets_msec_wd;
    dur_msec_sn=offsets_msec_sn-onsets_msec_sn;
    dur_msec_pr=offsets_msec_pr-onsets_msec_pr;
    
    
    save([expdir exp '/sound/onsets.mat'],'onsets_tr_wd','onsets_tr_sn','onsets_tr_pr','offsets_tr_wd','offsets_tr_sn','offsets_tr_pr','onsetsVectorMsec_wd','onsetsVectorMsec_sn','onsetsVectorMsec_pr',...
        'onsets_msec_wd','onsets_msec_sn','onsets_msec_pr','offsets_msec_wd','offsets_msec_sn','offsets_msec_pr','onsetsVectorTr_wd','onsetsVectorTr_sn','onsetsVectorTr_pr','dur_tr_wd','dur_tr_sn','dur_tr_pr','dur_msec_wd','dur_msec_sn','dur_msec_pr')
    % cols=flipud(gray(2));
    % imagesc([]); colormap gray;
end
%
% figure;
% for wdi=1:length(onsets_msec_wd);
%     line([onsets_msec_wd(wdi) onsets_msec_wd(wdi)]*1.5,[0 1],'color','k')
% end
% for sni=1:length(onsets_msec_sn);
%     line([onsets_msec_sn(sni) onsets_msec_sn(sni)]*1.5,[1 2],'color','k')
% end
% for pri=1:length(onsets_msec_pr);
%     line([onsets_msec_pr(pri) onsets_msec_pr(pri)]*1.5,[2 3],'color','k')
% end
% xlabel('time (sec)');
% set(gca,'fontsize',14,'ytick',[0.5:1:3],'yticklabels',{'Paragraph','Sentence','Word'});