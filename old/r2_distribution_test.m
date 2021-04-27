for i=1:1000;
x=rand(100,1); y1=rand(100,1);  y2=rand(100,1);
y1=y1-mean(y1);
y2=y2-mean(y2);
[b,~,~,~,stats]=regress(y1,x1);
r2_train(i,:)=stats(1);
sst=sum(y2.^2);
ssr=sum((y2-x*b).^2)
r2_test(i,1)=1-(ssr/sst);
end