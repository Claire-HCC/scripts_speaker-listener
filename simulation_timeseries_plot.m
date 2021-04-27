clear all
% loc='cluster';
close all
set_parameters;
timeUnit='tr' ;
%
%  froidir='restFc_isc15PercMasked_50Overlap_cluster6';
%  networks={'AUD','vLAN','dLAN','DMNa','DMNb'};

cols=jet(7);
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];

simN=30;
simi=1;
crop_start=5;
crop_end=0;
NRFs={'linear', 'log', 'triangle','time'};
gray=[0.7 0.7 0.7];
tr=1.5;


% filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr0.50_linear.mat',expdir,simN)));
filelists{1} = {'pause0.10_pauseL3.0_m3_v0.50_sr0.50_linear.mat'}'

for parami=1:length(filelists);
    filelist=filelists{parami};
    
    for fi=1:length(filelist);
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,filelist{fi}));
       
        
        figure('unit','centimeter','position',[0 0 15 6.5]);
        gdata=v_hrf_tr(:,:,simi);
        [levn,tn,~]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        
        gdata=gdata(:,keptT);
        gdata=zscore(gdata,0,2);
        for ni=1:levn;
            plot(keptT*1.5,gdata(ni,:),'color',cols(ni,:)*0.95,'linewidth',3)
            hold on;
        end
        % ylim([-2 5]);
        xlim([min(keptT) max(keptT)]*1.5);
        ylim([-4 3])
        %  line([0 0 ],get(gca,'ylim'),'color',gray);
        %  line(get(gca,'xlim'),[0 0],'color',gray);
        set(gca,'xtick',50:50:200,'ytick',[]);
        xlabel('Time (sec)');
        ylabel({'Simulate','BOLD response'})
        set(gca,'fontsize',14);
        % grid on
        hold off
    end
end
ylm=get(gca,'ylim');
rectangle('Position',[99 ylm(1) 10 ylm(2)-ylm(1)] ,'edgecolor','k','linewidth',3);
rectangle('Position',[10 ylm(1) 10 ylm(2)-ylm(1)] ,'edgecolor','k','linewidth',3);

%% zoom in
 set(gcf,'position',[0 0 5 10]);
title('')
ylabel('')
xlabel('')
set(gca,'xtick',[],'ytick',[])

% legend(networks,'orientation','horizontal','textcolor','w');
% legend boxoff;


