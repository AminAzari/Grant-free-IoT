%IHN
function Sc=decF(FD,D,y_pam_temp,ip,Fs,Tb,N,k,Nsamp,ncp,nzc,fc,com,FDTot,DTot,y0,cy0,IU,W,Nt,Np,zar,Th,TP)
t=[1/Fs:1/Fs:Tb];
if(length(D)==1)
    y_pam_temp_iij=y_pam_temp(floor(D*length(t)/Tb)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+floor(D*length(t)/Tb));
    pwO=zar^2;
elseif(strcmp(com,'sc'))
    for iu=1:length(D)
        Au(iu,:)=y_pam_temp(floor(D(iu)*length(t)/Tb)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+floor(D(iu)*length(t)/Tb));
    end
     
    [y_pam_temp_iij,cor]=scF(Au,FDTot,DTot,IU,W,Tb,Nsamp,ncp,nzc,N,k,Nt,Np,t,TP);
    pwO=sum(cor.*(zar'.^2)); 
elseif(strcmp(com,'mrc'))
    
    for iu=1:length(D)
        Au(iu,:)=y_pam_temp(floor(D(iu)*length(t)/Tb)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+floor(D(iu)*length(t)/Tb));
    end
    y_pam_temp_iij=mrcF(Au,FD,y0,cy0,fc,Fs,Tb,ncp,nzc,Nsamp,Nt);
elseif(strcmp(com,'egc'))
    
    for iu=1:length(D)
        Au(iu,:)=y_pam_temp(floor(D(iu)*length(t)/Tb)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+floor(D(iu)*length(t)/Tb));
    end
    y_pam_temp_iij=egcF(Au,FD,y0,cy0,fc,Fs,Tb,ncp,nzc,Nsamp,Nt);
    pwO=sum(zar.^2);
end
 
 No= (10^(-20.4)*200); 

pwSIN=sum(y_pam_temp_iij*y_pam_temp_iij')/length(y_pam_temp_iij); 
if(pwO/(abs(pwSIN+No-pwO))>Th) 
    Sc=1; 
else
    Sc=0;
end
% 
% estBit_noncoh= FuncDecP(t,y_pam_temp_iij,fc,Fs,Tb,FD,ncp,nzc,Nsamp);
% %--counting errors
% ipb=ip>mean(ip);
% nErr_pam_noncoh = size(find([ipb - estBit_noncoh]),2);  %counting the number of errors
% nErr=nErr_pam_noncoh/length(ip);
