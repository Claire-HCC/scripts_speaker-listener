clear all
% close all
loc='mypc';
set_parameters

iters=1000;
rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='merlin';
textSim=repmat([1 -1; -1 1],11,11);


for ri=2%1:length(rnames);
    rname=rnames{ri};
    load([expdir exp '/fmri/mat/roi/segment/segment_' role '_' rname '.mat']);
    
 data=zscore(data,0,2);
    for s=1:18;
        train_subjects=1:18;
        train_subjects=train_subjects(train_subjects~=s);
        train_data=mean(data(:,:,train_subjects),3);
        test_data=data(:,:,s);
        sim_temp=corr(train_data,test_data);
        
        sim(:,:,s)=sim_temp;
         ms(s,:)=grpstats(sim_temp(eye(22)==0), textSim(eye(22)==0),'mean');
    end
   
    acc=mean((ms(:,2)-ms(:,1))>0);
    
    for iter = 1:iters;
        shuff1 = randperm(22);
        shuff2 = randperm(22);
        
        for s=1:18;
            train_subjects=1:18;
            train_subjects=train_subjects(train_subjects~=s);
            train_data=mean(data(:,:,train_subjects),3);
            test_data=data(:,:,s);
            
            sim_temp=corr(train_data,test_data);
            sim_temp= sim_temp(shuff1, shuff1);
            shuffled_mriSim(:,:,s)=sim_temp;
            
           ms_null(s,:,iter)=grpstats(sim_temp(eye(22)==0), textSim(eye(22)==0),'mean');
            % Correlate
            
        end
    end
acc_null=mean(ms_null(:,2,:)-ms_null(:,1,:)>0,3);
    
    % calculate p-values
    p_acc= sum(acc_null>=acc)/iters;
    
    
    disp(sprintf('%s Checkerboard: acc=%0.3f (p=%0.3f)',rname,acc,p_acc));
end
%  clear r_pear r_spear r_pear_null r_spear_null p_pear p_spear p_thr_pear p_thr_spear

