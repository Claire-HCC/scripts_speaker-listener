clear all
close all
% loc='cluster';
% delete pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear.mat
% simulation_pallier(3,vars,speechRs,NRFs,pauseLens,pauseEffects);
% simulation_timeseries_plots

set_parameters;
timeUnit='tr' ;

cols=jet(7);
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];

simN=1;
crop_start=5;
crop_end=0;
NRFs={'linear', 'log', 'triangle','time'};
gray=[0.7 0.7 0.7];
tr=1.5;
simi=1;

xlm=[1 110];

% filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear.mat',expdir,simN)));
filelists{1} = {'pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear.mat'}'
fs1=[];
fs2=[];
for parami=1:length(filelists);
    filelist=filelists{parami};
    
    for fi=1:length(filelist);
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,filelist{fi}));
        [levn,tn,simN]=size(v);
        
        %%        plot simulated narrative structure
        gdata=onsetsVectors(:,:,simi);
        [levn,tn,~]=size(gdata);
        % keptT=(crop_start*1000*1.5+1):(tn-crop_end*1000*1.5);
        keptT=1:tn;
        figure('unit','centimeter','position',[0 0 23 9]);
        for ni=1:levn;
            plot(keptT/1000,(gdata(ni,keptT)==0)+(levn-ni+1),'color',cols(ni,:),'linewidth',ni*(4/6))
            hold on
        end
        xlim([min(keptT) max(keptT)]/1000);
        ylim([1 7.05]);
        set(gca,'fontsize',16,'ytick',[1:6]+0.5,'yticklabels', cellfun(@(x) ['Level ' x],cellstr(num2str(fliplr([1:6])')),'Uniformoutput',0),'xtick',[50 100],'xticklabels',[50 100]);
        xlabel('Time (sec)');
      title({'Simulated nested story structure:','consituent onsets'})
        
        %%        plot simulated neural activity
        gdata=v(:,:,simi);
        [levn,tn,~]=size(gdata);
        % keptT=(crop_start*1000*1.5+1):(tn-crop_end*1000*1.5);
        keptT=1:tn;
        for ni=1:levn;
            figure('unit','centimeter','position',[0 0 15 4]);
            plot(keptT/1000,gdata(ni,keptT),'color',cols(ni,:),'linewidth',ni*(4/6))
            axis tight
            xlim(xlim);
            ylim([min(gdata(:))-0.5 max(gdata(:))+0.5]);
            set(gca,'fontsize',12,'ytick',[],'xtick',[50 100],'xticklabels','');
            
            % axis off
            print(gcf,['temp' num2str(ni) '.tif'],'-dtiff','-r1000');
            fs1=[fs1; imread(['temp' num2str(ni) '.tif'])];
        end
        
        imshow(fs1);
        imwrite(fs1,'temp.tif')
        
        %% plot simulated BOLD response
        gdata=v_hrf_tr(:,:,simi);
        [levn,tn,~]=size(gdata);
        %    keptT=(crop_start+1):(tn-crop_end);
        keptT=(crop_start+1):(tn-crop_end);
        gdata=gdata(:,keptT);
        gdata=zscore(gdata,0,2);
        % gdata=zscore(gdata,0,2);
        for ni=1:levn;
            figure('unit','centimeter','position',[0 0 15 4]);
            
            plot(keptT*1.5,gdata(ni,:),'color',cols(ni,:)*0.95,'linewidth',ni*(4/6))
            
            % ylim([min(gdata(:))-1 max(gdata(:))+1]);
            xlim([min(keptT) max(keptT)]*1.5);
            ylim([-3 4])
            
            set(gca,'fontsize',12,'ytick',[],'xtick',[10 20 99 109],'xticklabels','');
       
            print(gcf,['tempp' num2str(ni) '.tif'],'-dtiff','-r1000');
            fs2=[fs2; imread(['tempp' num2str(ni) '.tif'])];
            
        end
        imshow(fs2);
        imwrite(fs2,'tempp.tif')
        
    end
end
% legend(networks,'orientation','horizontal','textcolor','w');
% legend boxoff;


fsize=[3 11.13];
figure('unit','centimeter','position',[ 0 0 fsize])
for ni=1:6; subplot(6,1,ni); text(0,0,['Level ' num2str(ni)],'fontsize',13); axis off ; end


