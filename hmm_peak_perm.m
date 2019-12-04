tic

clear all
close all
loc='cluster';
set_parameters;

froidir='mor';
rname='ABCshen222';

sample_gap=2:10;
permN=10;
minK=10;
for ei = 4;
    exp=experiments{ei};
    load([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm.mat'],'segments');
    
    % too few segments results in zeros data points for between-segment
    % correlation
    segments(1:(minK-1),:,:)=NaN;
    
    load([expdir '/' exp '/' 'fmri/timeseries/tr/roi/' froidir '/zscore_listenerAll_' rname '.mat'],'gdata');
    tn=size(gdata,2);
    
    % spatial correlation be- tween all pairs of time points that were separated by four time points
    C_sample=zeros(tn,tn);
    
    for ti=1:(tn-sample_gap);
        C_sample(ti,min(ti+sample_gap,tn))=1;
    end
    
    for s=1:size(gdata,3);
        self=gdata(:,:,s);
        C =corr(self);
        % made symmetric
        C(eye(size(C,1))==1)=0;                        % put 1 on the diagonal
        Cs(:,s)=C(:);
    end
    
    % Within- v.s. Between segment pattern Correlation
    wbr=nan(size(segments,1),size(gdata,3));
    wbr_perm=nan(size(segments,1),size(gdata,3),permN);
    for Ki=minK:size(segments,1);
        evi=segments(Ki,:,s);
        
        C_withini=(repmat(evi,tn,1)-repmat(evi',1,tn))==0;
        C_betweeni=(C_withini==0);
        % within vs. cross boundary pattern correlation, following Baldassano (2017)
        Cs_within=Cs(C_withini(:)==1 & C_sample(:)==1,:);
        Cs_between=Cs(C_withini(:)~=1 & C_sample(:)==1,:);
        
        if size(Cs_within,1)>5 & size(Cs_between,1)>5;
            wbr(Ki,:)= nanmean(Cs_within,1)-nanmean(Cs_between,1);
        else
            wbr(Ki,:)=NaN;
        end
    end
    
    % perm
    for Ki=minK:size(segments,1);
        if sum(isnan( wbr(Ki,:)))==0;
            for pi=1:permN;
                evi_perm=zeros(size(evi));
                evi_perm(randsample(length(evi),max(evi)-1))=1;
                evi_perm=cumsum(evi_perm)+1;
                C_withini_perm=(repmat(evi_perm,tn,1)-repmat(evi_perm',1,tn))==0;
                C_betweeni_perm=(C_withini_perm==0);
                % within vs. cross boundary pattern correlation, following Baldassano (2017)
                Cs_within_perm=Cs(C_withini_perm(:)==1 & C_sample(:)==1,:);
                Cs_between_perm=Cs(C_withini(:)~=1 & C_sample(:)==1,:);
                Cs_between_perm=Cs(C_withini_perm(:)~=1 & C_sample(:)==1,:);
                wbr(Ki,:)= nanmean(Cs_within_perm,1)-nanmean(Cs_between_perm,1);
            end
        end
    end
end

save([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm.mat'],'segments','wbr','wbr_perm');
toc