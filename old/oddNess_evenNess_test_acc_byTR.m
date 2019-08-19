clear all
% close all
loc='mypc';
set_parameters

iters=100;
rnames={'HG_L','vPCUN','STC_L'};
role='listener';
exp='merlin';
timeUnit='tr';

load([expdir exp '\sound\listenerEvents.mat']);
textV=-ones(1,length(listenerEvents_trVector));
textV(mod(listenerEvents_trVector,2)==1)=1;
textV=textV(listenerEvents_trVector>0);
textSim=(textV')*(textV);
eventn=size(textSim,2);

for ri=1:3;%:length(rnames);
    rname=rnames{ri};
    load([expdir exp '/fmri/mat/roi/' timeUnit '/' timeUnit '_' role '_' rname '.mat']);
    
    % delete pause scane
    data=zscore(data,0,2);
    data=data(:,listenerEvents_trVector>0,:);
    
    
    for s=1:18;
        train_subjects=1:18;
        train_subjects=train_subjects(train_subjects~=s);
        train_data=mean(data(:,:,train_subjects),3);
        test_data=data(:,:,s);
        sim_temp=corr(train_data,test_data);
        sim(:,:,s)=sim_temp;
        segn=size(sim_temp,1);
        ms(s,:)=grpstats(sim_temp(eye(segn)==0), textSim(eye(segn)==0),'mean');
    end
    acc=mean((ms(:,2)-ms(:,1))>0);
    
    for iter = 1:iters;
        shuff1 = randperm(eventn);
        shuff2 = randperm(eventn);
        for s=1:18;
            train_subjects=1:18;
            train_subjects=train_subjects(train_subjects~=s);
            train_data=mean(data(:,:,train_subjects),3);
            test_data=data(:,:,s);
            sim_temp=corr(train_data,test_data);
            sim_temp=sim_temp(shuff1, shuff1);
            shuffled_mriSim(:,:,s) = sim_temp;
            ms_null(s,:,iter)=grpstats(sim_temp(eye(segn)==0), textSim(eye(segn)==0),'mean');
            % Correlate
            
        end
    end
    acc_null=mean(ms_null(:,2,:)-ms_null(:,1,:)>0,3);
    
    p_acc = sum( acc_null>=acc)/iters;
    
    disp(sprintf('%s Checkerboard: acc=%0.3f (p=%0.3f)',rname,acc,p_acc));
end
%  clear r_pear r_spear r_pear_null r_spear_null p_pear p_spear p_thr_pear p_thr_spear

