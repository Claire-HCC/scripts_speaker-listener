timeUnit='tr';
froidir='mor';
lags=-10:10;
for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','peak','peakLags','peakLags_pfdr');
    peakLags_pfdr_all(:,ei)=peakLags_pfdr;
    peakLags_all(:,ei)=peakLags;
    
end
% 
% 
% for ei1=1:4;
%     for ei2=1:4;
%         ris=find( ~isnan(peakLags_pfdr_all(:,ei1)) & ~isnan(peakLags_pfdr_all(:,ei2)));
%         [r(ei1,ei2) p(ei1,ei2)]=corr(peakLags_pfdr_all(ris,ei1),peakLags_pfdr_all(ris,ei2),'type','spearman');
%         n(ei1,ei2)=length(ris);
%     end
% end

ris=find(sum(isnan(peakLags_all),2)==0);
[r p]=corr(peakLags_all(ris,:),'type','spearman');

% r matrix
figure;
imagesc(r,[-1 1]); 
set(gca,'xtick',1:4,'xticklabels',experiments)
set(gca,'ytick',1:4,'yticklabels',experiments)
colorbar
title({'Peak lags', 'Spearman correlation'});
for ei1=1:4;
    for ei2=1:4;
        if p(ei1,ei2)<.05; star='*'; else star=''; end;
        text(ei1,ei2,[num2str(round(r(ei1,ei2),2)) star],'fontsize',16,'HorizontalAlignment','center');
    end
end
set(gca,'fontsize',16);

% scatter plot
figure; 
ei1=3; ei2=4;
text(peakLags_all(:,ei1)+rand(61,1)-0.5,peakLags_all(:,ei2)+rand(61,1)-0.5,rnames,'HorizontalAlignment','center');
title(sprintf('R=%.2f; p=%.3f',r(ei1,ei2),p(ei1,ei2)));
xlim([-11 11]);
ylim([-10 10]);
xlabel(experiments{ei1});
ylabel(experiments{ei2});
set(gca,'fontsize',16);
grid on