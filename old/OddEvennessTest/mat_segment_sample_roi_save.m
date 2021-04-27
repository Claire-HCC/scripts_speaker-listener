clear all
set_parameters;
sampleN=2;% (merlin: mean 25 scans (range 5-92); sherlock mean 30 (range 4-76))

rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='sherlock';

load([expdir exp '\sound\listenerEvents.mat']);
labels=listenerEvents_trVector;

for ri=1:length(rnames);
    rname=rnames{ri};
    load([expdir exp '/fmri/mat/roi/tr/tr_' role '_' rname '.mat']);
    data=zscore(data,0,2);
    
    
    samplei=mod(1:length(labels),sampleN)';
    n=1;
    for segi=1:max(labels);
        for sample=1:sampleN;
            i=(labels==segi & samplei==(sample-1));
            data2(:,n,:)=mean(data(:,i,:),2);
            labels2(n)=mean(labels(i));
            n=n+1;
        end
    end
    
    gdata=data2;
    labels=labels2';
    fname = sprintf('%s\\%s\\fmri\\mat\\roi\\segment\\segment_%s_sample%d_%s.mat',expdir,exp,role,sampleN,rname);
    % save(fname,'gdata','labels')
    
    % figrue
    segn=length(labels2);
    dat=mean(data2,3);
    d=pdist(dat','euclidean');
    Y=mdscale(d,2);
    
    color=zeros(segn,3);
    color(mod(labels2,2)==1,:)=repmat([1 0 0],sum(mod(labels2,2)==1),1);
    color(mod(labels2,2)==0,:)=repmat([0 0 1],sum(mod(labels2,2)==0),1);
    color=color.*repmat([0.3:(0.7/segn):(1-0.7/segn)]',1,3);
    
    figure;
    scatter(Y(:,1),Y(:,2),100,color,'filled');
    text(Y(:,1),Y(:,2), cellstr(num2str(labels2')),'horizontal','left', 'vertical','bottom','fontsize',10);
    title({rname, 'Group average' , 'sample per segment:' num2str(sampleN)});
    
    print(gcf,[expdir '/graph/OddEveness_test/' exp '_segment_sample' num2str(sampleN) '_' rname '.png'],'-dpng')
    
    clear data data2
end
