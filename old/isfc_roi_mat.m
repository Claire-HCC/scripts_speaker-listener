% role1 roi to role2 wholebrain
% negative peak T means role1 precedes, positive peakT means role2 precedes
loc='mypc';
set_parameters
%
role1='listener01';
role2='listener01_others';

% role1='speaker01';
% role2='listenerZscoreMean';%'listener01_others';

relation='LL_s01_leave1out';%'SL';
lags=[-10:10]; %scan

p_thr=0.05;

iters=1000;
roi_selected={'PMC_L','HG_L','STC_L','pIFG_L','aANG_L','vPCUN','vmPFC'};
roiN=roiN;
for ei=1%:2;%1:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat/roi/' exp '_' role1 '_rois.mat' ],'data');
    data1=table2array(data(:,roi_selected));
    
    load([expdir experiments{ei} '/fmri/mat/roi/' exp '_' role2 '_rois.mat' ],'data');
    data2=table2array(data(:,roi_selected));
    
    for i=1:roiN;
        for j=1:roiN;
            lagcc_temp= lagcorr(data1(:,i),data2(:,j),lags)';
            [peakR peakT]=max((lagcc_temp)');
            peakRmat(i,j)=peakR;
            peakTmat(i,j)=lags(peakT);
            
        end
    end
    
    save([expdir exp '/fmri/mat/roi/' role1 '_' role2 '.mat' ],'peakRmat','peakTmat','roi_selected','-v7.3');
    
    for iter=1:iters;
        data1_perm=phase_rand(data1,1);
        data2_perm=phase_rand(data2,1);
        r=corr(data1_perm,data2_perm);
        rmat_perm(:,:,iter)=r;
    end
    save([expdir exp '/fmri/mat/roi/' role1 '_' role2 '.mat' ],'peakRmat','peakTmat','rmat_perm','roi_selected','-v7.3');
end


figure;
% peakTmat(tril(ones(size(peakTmat)),-1)==1)=NaN;
p=sum((peakRmat<rmat_perm),3)/(iters);
sig_fwe= p < (p_thr/(roiN^2*length(lags)));
sig_fdr=reshape(fdr0(p(:),p_thr),roiN,roiN);
peakTmat_thresholded=peakTmat;
peakTmat_thresholded(sig_fdr==0)=NaN;

imagesc(repmat(-20,roiN,roiN));
hold on
m=imagesc(peakTmat_thresholded,[-10 10]);
set(m,'AlphaData',~isnan(peakTmat_thresholded));
hold off
colormap('hot')
colorbar
set(gca,'xtick',1:roiN);
set(gca,'xticklabel',roi_selected);
set(gca,'ytick',1:roiN);
set(gca,'yticklabel',roi_selected);
xlabel(role2);
ylabel(role1);
%  colorbar
% black if R value lower than threshold

[x, y] = meshgrid(1:roiN);  % Create x and y coordinates for the strings
 hStrings = text(x(:), y(:), num2str(peakTmat_thresholded(:)),'HorizontalAlignment', 'center','fontsize',12);





