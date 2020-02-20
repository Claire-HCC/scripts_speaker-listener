clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:-4};
binSize=10;

for ei=4;%1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for si=1:48;
            if exist(([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat']));
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat'],'r2');
                keptT=[11:(size(r2,2)-12)];
                r2_perm(:,:,si)=r2;
                
                for ri=1:length(rnames);
                    b=polyfit(keptT',r2_perm(ri,keptT,si)',1);
                    pv_perm(ri,:,si) = polyval(b,keptT);
                end
                
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag-3-3_permSL' num2str(si) '.mat'],'r2');
                r2_perm0(:,:,si)=r2;
                
                for ri=1:length(rnames);
                    b=polyfit(keptT',r2_perm0(ri,keptT,si)',1);
                    pv0_perm(ri,:,si) = polyval(b,keptT);
                end
            else
                pv0_perm(:,:,si)=NaN;
                pv_perm(:,:,si)=NaN;
            end
        end
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag-3-3'  ],'rnames','r2');
        r20=r2;
        for ri=1:length(rnames);
            b=polyfit(keptT',r20(ri,keptT)',1);
            pv0(ri,:) = polyval(b,keptT);
        end
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'rnames','r2');
        
        for ri=1:length(rnames);
            b=polyfit(keptT',r2(ri,keptT)',1);
            pv(ri,:) = polyval(b,keptT);
        end
        
    end
end


for ri=25;%1:61;
    figure;
subplot(1,2,1);
ciplot_claire(squeeze(pv_perm(ri,:,:))',keptT,'k',0.3)
hold on;
plot(keptT,pv(ri,:),'r')
hold off
grid on;
ylim([0 0.5])
title({exp,rnames{ri},'lag-10~-4 TR'})

subplot(1,2,2);
ciplot_claire(squeeze(pv0_perm(ri,:,:))',keptT,'k',0.3)
hold on;
plot(keptT,pv0(ri,:),'r')
hold off
grid on;
ylim([0 0.5])
title({exp,rnames{ri},'lag-3~3 TR'})
pause
close gcf
end
