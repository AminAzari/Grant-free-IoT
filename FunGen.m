function [y,ip,y0,cy0,zar]=FunGen(N, Nsamp,Fs,Tb,fc,fd,d0,ncp,nzc,Nt,k,cons,dis,Pt,ph )

%generating random binary sequence
ip = cons+(rand(1,N)>0.5);  %generating 0,1 with equal probability
pam_symbols_temp1 = ip(1:N/2);
pam_symbols_temp2 = ip(N/2+1:N);
pam_symbols_temp=pam_symbols_temp1+1i*pam_symbols_temp2;
%scatter(real(dqpsk_symbols), imag(dqpsk_symbols))
%white gaussian noise, 0 mean
pam_symbols_te = rectpulse(pam_symbols_temp,Nsamp);

Ase= transpose(lteZadoffChuSeq(1,nzc));
Ase= rectpulse(Ase,Nsamp);
pam_symbols_t = [pam_symbols_te(end-ncp+1:end) Ase    pam_symbols_te ];
%dqpsk_symbols = [seq dqpsk_symbols_t ];

%%----------------TRansmitter----   Modulation to carrier freq
t=[1/Fs:1/Fs:Tb];
y0=[];
for i=1:length(pam_symbols_t)
    y0=[y0 (1*pam_symbols_t(i)*exp(1i*2*pi*(fc+fd)*[(i-1)*Tb+1/Fs:1/Fs:i*Tb]))];
end
py=sum(y0*y0')/length(y0);
y0=sqrt(1/py)*y0; 
% py=sum(y0*y0')/length(y0);
y=zeros(1, Nt*length(t)*1/Tb); 
%+1i*randn(1,Nt*length(t));


for iu=1:length(d0)
    PL=133 + 38.3*log(dis./1000)+10*log10(exprnd(1.775,1,1));
    y01 = sqrt((10^(Pt/(10)))/(10^(PL/10)) )*y0;
    zar(iu)=sqrt((10^(Pt/(10)))/(10^(PL/10)) );
    
% dis
    %     [length(floor(d0(iu)*length(t))+1:floor(length(t)*(N/k*Nsamp+nzc*Nsamp+ncp+d0(iu)))) length(y0)]
    
    y(floor(d0(iu)*length(t)/Tb)+1:length(t)*(N/k*Nsamp+nzc*Nsamp+ncp)+floor((length(t)*d0(iu)/Tb)) )=y01;
end
 
cy0=fd;
