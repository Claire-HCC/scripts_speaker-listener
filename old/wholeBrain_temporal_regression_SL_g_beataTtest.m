
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

connames={'StephensSprecedes','StephensSync','StephensLprecede',...
    'mySprecedes','myLprecede'};

% lags=-10:10;
% contrasts={[12/9*ones(1,9) -ones(1,3) -ones(1,9)],...
%     [-ones(1,9) 18/3*ones(1,3) -ones(1,9)],...
%     [-ones(1,9) -ones(1,3) 12/9*ones(1,9)],...
%     [ ones(1,9) zeros(1,3) -ones(1,9) ],...
%     [-ones(1,9) zeros(1,3) ones(1,9) ]};

lags=-4:4;
contrasts={[2*ones(1,3) -ones(1,3) -ones(1,3)],...
    [-ones(1,3) 2*ones(1,3) -ones(1,3)],...
    [-ones(1,3) -ones(1,3) 2*ones(1,3)],...
    [ ones(1,3) zeros(1,3) -ones(1,3) ],...
    [-ones(1,3) zeros(1,3) ones(1,3) ]};

for ei=1%:4;
    exp=experiments{ei};
   % load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'pfdr');
    load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','keptvox');
    [~,~,listenerN]=size(b);
    
    for ci=1:length(contrasts);
        contrast=contrasts{ci};
        conname=connames{ci};
        c=nan([length(keptvox) listenerN]);
        
        subjs=find(~isnan(squeeze(b(1,1,:))));
        for si=1:length(subjs);
            s=subjs(si);
            c(:,s)=b(:,:,s)*contrast';
        end
        [~,p,~,stats]=ttest(c',0,'tail','right');
        t=stats.tstat;
        [~,~,pfdr]=fdr(p);
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' conname ],'c','p','t','keptvox','b','pfdr');
        mat=nan([voxn 1]);
        
        mat(keptvox(pfdr<.05))=t(pfdr<.05);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' conname '_t_pfdr.nii']);
    end
end

%% plot
cols={'b','y','r','b','r'};
lags_tested={-4:4,-10:10};
fsize=[10 45];
for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    figure('unit','centimeter','position',[0 0 fsize]);
    for ei=1:4;
        exp=experiments{ei};
        subplot(1,4,ei)
        
        for ci=1:3;
            contrast=contrasts{ci};
            conname=connames{ci};
            
            load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' conname ],'c','p','t','keptvox','b','pfdr');
            plot(lags,squeeze(nanmean(nanmean(b(pfdr<.05,:,:),3))),'linewidth',2,'color',cols{ci});
            % ciplot_claire(squeeze(nanmean(b(pfdr<.05,:,:),1))',lags,cols{ci},0.3);
            n{ci}=sum(pfdr<.05);
            set(gca,'fontsize',14);
            hold on
        end
        grid on;
        title(exp);
        % legend(cellfun(@(x,y) sprintf('%s (voxel N=%d)',x,y),connames(1:3),n(1:3),'UniformOutput',0));
        legend(cellfun(@(x) sprintf('voxel N=%d',x),n(1:3),'UniformOutput',0));
        legend boxoff
        hold off
    end
    xlabel('Lags (TR)');
    ylabel('Averaged beta');
    
    figure('unit','centimeter','position',[0 0 fsize]);
    for ei=1:4;
        exp=experiments{ei};
        subplot(1,4,ei)
        for ci=[4 5]
            contrast=contrasts{ci};
            conname=connames{ci};
            
            load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' conname ],'c','p','t','keptvox','b','pfdr');
            plot(lags,squeeze(nanmean(nanmean(b(pfdr<.05,:,:),3))),'linewidth',2,'color',cols{ci});
            % ciplot_claire(squeeze(nanmean(b(pfdr<.05,:,:),1))',lags,cols{ci},0.3);
            n{ci}=sum(pfdr<.05);
            set(gca,'fontsize',14);
            hold on
        end
        grid on;
        title(exp);
        %   legend(cellfun(@(x,y) sprintf('%s (voxel N=%d)',x,y),connames(4:5),n(4:5),'UniformOutput',0));
        legend(cellfun(@(x) sprintf('voxel N=%d',x),n(4:5),'UniformOutput',0));
        legend boxoff
        hold off
    end
    xlabel('Lags (TR)');
    ylabel('Averaged beta');
end