clear all
close all

loc='cluster';
set_parameters;
exp='sherlock';
role='listener';
sampleN=2;
dmethod='euclidean';
iters=1000;
p_thr=0.05;
rnames={'HG_L','vPCUN','STC_L'};

for ri=1:length(rnames);
    
    set_parameters;
    rname=rnames{ri};
    
    fname = sprintf('%s/%s/fmri/mat/roi/segment/segment_%s_sample%d_%s.mat',expdir,exp,role,sampleN,rname);
    load(fname);
    gdata=zscore(gdata,0,2);
    
    labels_abc=zeros(length(labels),1);
    labels_abc(find(mod(labels,2)==1 & labels>0 ))=1;
    labels_abc(find(mod(labels,2)==0 & labels >0))=2;
    
    segn=sum(labels>0);
    seg_inds=find(labels>0);
    
    cbs=nchoosek(1:size(gdata,3),2); % test set
    cbn=size(cbs,1);
    
    for cbi=1:size(cbs,1);
        subjects=1:size(gdata,3);
        subjs_test=cbs(cbi,:);
        subjs_train=subjects(~ismember(subjects,cbs(cbi,:)));
        
        subjSet_train= mean(gdata(:,:,subjs_train),3);
        subjSet_test= mean(gdata(:,:,subjs_test),3);
        
        svm_model = fitcsvm(subjSet_train(:,seg_inds,:)',labels_abc(seg_inds));
        [label_predicted(1,:,cbi)]=predict(svm_model ,subjSet_test');
    end
    
    
    for iter=1:iters;
        
        segment_shuffled=randperm(segn);
        
        for cbi=1:size(cbs,1);
            
            subjects=1:size(gdata,3);
            subjs_test=cbs(cbi,:);
            subjs_train=subjects(~ismember(subjects,cbs(cbi,:)));
            
            subjSet_train= mean(gdata(:,seg_inds(segment_shuffled),subjs_train),3);
            subjSet_test= mean(gdata(:,:,subjs_test),3);
            
            svm_model = fitcsvm(subjSet_train',labels_abc(seg_inds));
            label_predicted_perm(1,:,cbi,iter)=predict(svm_model ,subjSet_test');
            
        end
    end
    
    
    %% acc graph
    clear labels_temp;
    labels_temp(1,:,:)=repmat(labels_abc(seg_inds),1,cbn);
    acc=(label_predicted(:,seg_inds,:)==labels_temp);
    acc_seg=grpstats(squeeze(acc),labels(seg_inds),'mean');
    m=squeeze(mean(acc,3));
    m_seg=grpstats(m,labels(seg_inds),'mean');
    ci95=ci(acc_seg',0.95)
    
    clear labels_temp;
    labels_temp(1,:,:,:)=repmat(labels_abc(seg_inds),1,cbn,iters);
    acc_perm=(label_predicted_perm(:,seg_inds,:,:)==labels_temp);
    m_perm=squeeze(mean(acc_perm(:,seg_inds,:,:),3));
    m_perm_seg=grpstats(m_perm,labels(seg_inds),'mean');
    p_pos=sum(m_perm_seg>repmat(m_seg,1,iters),2)/iters;
    p_neg=sum(m_perm_seg<repmat(m_seg,1,iters),2)/iters;
    sig=fdr0(p_pos,p_thr/2)+fdr0(p_neg,p_thr/2)
    
    figure
    subplot(1,4,1);
    bar(mean(m_seg),'k');
    hold on
    plot([0 2],[mean(m_perm_seg(:)) mean(m_perm_seg(:))],'--r');
    if sum(mean(m_perm_seg,1)>mean(m_seg))/iters<p_thr;
        text(1,mean(m_seg)+0.03,'*','horizontal','center','fontsize',18);
    end
    hold off
    xlim([0 2]);
    ylim([0 1]);
    set(gca,'xtick',[]);
    xlabel({'All segment'});
    ylabel({'Classification Accuracy'})
    title(rname);
    
    subplot(1,4,[2:4]);
    bar(1:length(m_seg),m_seg,'k');
    hold on
    % errorbar(m_seg,ci95);
    plot(1:length(m_seg),mean(m_perm_seg,2),'--r');
    
    text(find(sig==1),m_seg(sig==1)+0.01,'*','horizontal','center','fontsize',18);
    hold off
    set(gca,'xtick',1:length(m_seg),'fontsize',7);
    xlim([0 length(m_seg)+1]);
    ylim([0 1]);
    xlabel('Individual segment');
    fname=sprintf('%s/graph/fMRI_segment/roi_ABacc_sample%d_cbn%d_iters%d/%s.png',expdir,sampleN,cbn,iters,rname)
    print(gcf,[expdir '/graph/OddEveness_test/' exp '_segment_sample' num2str(sampleN) '__' rname '_classification.png'],'-dpng')
  %  clearvars -except rnames label iters dmethod p_thr sampleN exp role
end
beep

