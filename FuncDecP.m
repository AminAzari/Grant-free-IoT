%%---------------------data decoding
function estBit_noncoh=FuncDecP(t,y_pam_temp,fc,Fs,Tb,cfo,ncp,nzc,Nsamp)
y_rx=[];
y_r=[];
y_rk=[];
y_rj=[];
p_rk=[];
 
j=0;
for i=1:length(t):length(y_pam_temp)
    j=j+1;
    AAk=y_pam_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc+cfo)*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
    AAj=y_pam_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc)*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
    y_rk=[y_rk AAk ];
    y_rj=[y_rj AAj ];
    if(i<=(nzc*Nsamp+ncp)*length(t))
        AAx=y_pam_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc+cfo)*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
        y_rx=[y_rx AAx ];
    else
        AAr=y_pam_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc+cfo)*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
        y_r=[y_r mean(AAr) ]; 
    end
    
end


% figure()
% 
% FX=[1:1:Fs-1];
% periodogram(y_rj(1:Fs*Tb),[],FX,Fs);
% figure()
% periodogram(y_rj(1:Nsamp*Fs*Tb),[],FX,Fs);
% 
% figure()
% periodogram(y_rj,[],FX,Fs);



%     Ffil =firceqrip(40,2*75/Fs,[0.05 0.03]);
 %     y_pam = intdump(y_nj,Nsamp);
y_pam = intdump(y_r,Nsamp);
ry_r=real(y_pam);
ry_i=imag(y_pam);
eyr=ry_r>mean(ry_r);
eyi=ry_i>mean(ry_i);
estBit_noncoh=[eyr eyi];

