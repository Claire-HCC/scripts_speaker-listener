gf=fft(g); 
for fi=1:floor(length(g)/2);
gf_temp=zeros(size(gf));
gf_temp(end-fi)=gf(end-fi);gf_temp(fi+1)=gf(fi+1);
gf_if(fi,:)=(ifft(gf_temp,'symmetric'));
end;