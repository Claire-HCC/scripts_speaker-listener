set_parameters;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));
%  nii=roiTable2wholeBrainNii_mor([roi_ids, [1:61]]);
%              save_nii(nii,[expdir '/roi_mask/mor/roi_order.nii']);

rgb=[];

% yellolw
ris=find(strcmp('DMN1',roi_table.network));
rgb_temp=[];
hsv=[];
h=[50 60]/360;
s=[1 0.5];
v=[0.8 1];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end

rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);

rgb(ris,:)=rgb_temp(1:length(ris),:);

%%turquoise 
ris=find(strcmp('Visual',roi_table.network));
rgb_temp=[];
hsv=[];
h=[166 177]/360;
s=[1 0.6];
v=[0.6 1];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end
rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);

rgb(ris,:)=rgb_temp(1:length(ris),:);

% red 
rgb_temp=[];
hsv=[];
ris=find(strcmp('Auditory_Language',roi_table.network));
h=[0 5]/360;
s=[1 0.6];
v=[0.6 1];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end
rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);

rgb(ris,:)=rgb_temp(1:length(ris),:);

% orange
ris=find(strcmp('DMN2',roi_table.network));
rgb_temp=[];
hsv=[];
h=[20 30]/360;
s=[0.7 1];
v=[0.7 1];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end
rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);

rgb(ris,:)=rgb_temp(1:length(ris),:);

% blue
ris=find(strcmp('Attention',roi_table.network));
rgb_temp=[];
hsv=[];
h=[230 240]/360;
s=[0.9 0.6];
v=[1 0.7];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end
rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);

rgb(ris,:)=rgb_temp(1:length(ris),:);


% green
ris=find(strcmp('Executive',roi_table.network));
rgb_temp=[];
hsv=[];
h=[100 110]/360;
s=[1 0.5];
v=[0.7 1];
rn=length(ris);
for ri=1:rn;
hsv(end+1,:)=[h(1)+(ri-1)*(h(2)-h(1))/rn s(1)+(ri-1)*(s(2)-s(1))/rn v(1)+(ri-1)*(v(2)-v(1))/rn];
end
rgb_temp=hsv2rgb(hsv);
colormap(rgb_temp);
colorbar
rgb_temp=rgb_temp(randperm(size(rgb_temp,1)),:);
rgb(ris,:)=rgb_temp(1:length(ris),:);


colormap(rgb);
colorbar

roi_table.rgb=rgb;
save([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

% xjview
% colormap([gray(64); rgb; repmat([0 0 0],3,1)]);
% set max & min as 65 1
