% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
close all

loc='mypc';
set_parameters;
timeUnit='tr';
k=6;
cutoff=0.1;
alpha_cutoff=0.01;

eis=[1 2 4 11 12 9 10 13];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/wholeBrain_psd'  ],'Sxx','keptvox','keptT','freq','win');
    [idx, C]= kmeans(abs(Sxx(:,freq<cutoff)),k,'Replicates',20,'Distance','sqeuclidean');
    freq=freq(freq<cutoff);
    
    alpha=sum(C(:,freq<alpha_cutoff),2);
    [~,idx_order]=sort(alpha);
    idx_new=zeros(size(idx));
    for ki=1:k;
        idx_new(idx==ki)=idx_order(ki);
    end
    idx=idx_new;
    C=C(idx_order,:);
    alpha=alpha(idx_order);
    save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/wholeBrain_psd_cluster'  ],'C','idx','keptvox','keptT','freq','win','alpha');
    
    for ki=1:k;
        roimask=zeros(voxn,1);
        roimask(keptvox,1)=(idx==ki);
        nii=mat2nii(roimask);
        save_nii(nii,[expdir '/roi_mask/psd_clusters/nii/' exp '_cluster' num2str(ki) '.nii'  ])
        save([expdir '/roi_mask/psd_clusters/mat/' exp '_cluster' num2str(ki) '.mat'  ],'roimask')        
    end
end


cols=jet(7);
eis=[1 2 4 11 12 9 10 13];
fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
   load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/wholeBrain_psd_cluster'  ],'C','idx','keptvox','keptT','freq','win');
  
    subplot(3,3,eii)
    
    for ni=1:6;
        plot(log(freq),C(ni,:),'color',cols(ni+1,:),'linewidth',2);
        hold on
    end
    
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')])
    hold off
    xlabel('Frequency (Hz)');
    ylabel('PSD (1/Hz)');
    xlim(log([freq(2) freq(end)]))
    ylim([0 50])
    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
  
    grid on
    set(gca,'fontsize',14)
end


eis=[1 2 4 11 12 ];
mat=zeros(voxn,1);
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
   load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/wholeBrain_psd_cluster'  ],'C','idx','keptvox','keptT','freq','win');
   mat(keptvox,eii)=idx;
end

vis=(sum(mat~=0,2)==5)
mat2=mat(vis,:)
[~,sorti]=sort(mat2(:,3));
imagesc(mat2(sorti,:))