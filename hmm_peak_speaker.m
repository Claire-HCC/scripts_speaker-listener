tic

clear all

close all
% loc='cluster';
set_parameters;

froidir='mor';
rname='dPCC';
sample_gap=2:10;

fsize=[38 10];
figure('unit','centimeter','paperposition',[0 0 fsize],'position',[0 0 fsize],'papersize',fsize);
for ei = 3;
    exp=experiments{ei};
    load([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm_speaker.mat'],'segments');
    % too few segments results in zeros data points for between-segment
    % correlation
    segments(1:9,:,:)=NaN;
    
    load([expdir '/' exp '/' 'fmri\timeseries\tr\roi\' froidir '/zscore_speaker_' rname '.mat'],'data');
    tn=size(data,2);
    
    % spatial correlation be- tween all pairs of time points that were separated by four time points
    C_sample=zeros(tn,tn);
    
    for ti=1:(tn-sample_gap);
        C_sample(ti,min(ti+sample_gap,tn))=1;
    end
    
    % Within- v.s. Between segment pattern Correlation
    wbr=nan(size(segments,1),size(data,3));
    
    for s=1:size(data,3);
        self=data(:,:,s);
        C =corr(self);
           % made symmetric
        C(eye(size(C,1))==1)=0;                        % put 1 on the diagonal
        
        for Ki=1:size(segments,1);
            evi=segments(Ki,:,s);
            
            C_withini=(repmat(evi,tn,1)-repmat(evi',1,tn))==0;
            C_betweeni=(C_withini==0);
            % within vs. cross boundary pattern correlation, following Baldassano (2017)
            C_within=C(C_withini==1 & C_sample==1);
            C_between=C(C_withini~=1 & C_sample==1);
            
            if length(C_within)>5 & length(C_between)>5;
                wbr(Ki,s)= nanmean(C_within(:))-nanmean(C_between(:));
            else
                wbr(Ki,s)=NaN;
            end
        end
        
    end
    %                figure; histogram(C_between(:),'BinEdges',-1:0.1:1,'Normalization','probability');
    %            hold on;
    %            histogram(C_within(:),20,'facecolor','r','BinEdges',-1:0.1:1,'Normalization','probability');
    %            hold off
    %            title(['Segment N: ' num2str(Ki)]);
    %
    % segment number is limited by story length.
    % short story with large segment number results in little sample points
    % for within-segment correlation
    K= find(sum(isnan(wbr(:,sum(~isnan(wbr),1)>1)),2)<1);
   
    segmentLength =size(data,2)*tr(ei)./K;
    
   subplot(1,4,ei)
    %  subplot(1,2,1);
    wbr_narm=wbr(K,sum(~isnan(wbr),1)>1);
    %     imagesc(zscore(wbr_narm,0,1)');
    %
    %     set(gca,'xticklabels',K(  get(gca,'xtick')));
    %     xlabel('segment number');
    %     ylabel('subject')
    %     title({exp,rname});
    
    % subplot(1,2,2);
    
    ci95=ci(wbr_narm',0.95);
    m=nanmean(wbr_narm,2)';
    
    KperMin=K/(tr(ei)*size(data,2)/60);
    
    ciplot(m-ci95/2,m+ci95/2,KperMin,'k',0.3);
    hold on;
    
    plot(KperMin,m,'b','linewidth',1.5);
    hold off
    grid on
    
    % find segment number with the largest wbr among segment number with
    % wbr significantly larger than zero
    try
    [mx,mxi]=max(m-ci95/2);
    peakK=K(mxi);
    peakKpermin=KperMin(mxi);
    peak_segmentLength= size(data,2)*tr(ei)/peakK;
    catch
        peakK=NaN;
        peak_segmentLength=NaN;
    end
        title({exp,rname,sprintf('peak segment legnth: %.02f sec', peak_segmentLength) ,sprintf('peak segment number/min: %.2f', peakKpermin),['K=' num2str(peakK)]}); %
        
    xlabel('Segment Number / min');
    ylabel({'within- minus between-','boundary pattern correlation (r)'});
    xlim([min(KperMin) max(KperMin)]);
    line([min(KperMin) max(KperMin)],[0 0],'color','k','linestyle','--');
    
    %   imagesc(diff(squeeze((segments(peakK,:,:))),1)'); colormap gray
end
