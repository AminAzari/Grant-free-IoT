%IHN
function y =remF(zar,y,dSuc,D,FD,Fs,Tb,ipin,N,nzc,ncp,Nsamp,Np,fc,k,fas,H,ncz)

t=[1/Fs:1/Fs:Tb];
for ic=1:Np
    if(dSuc(ic)==1)
        ip=ipin(ic,:);
        pam_symbols_temp1 = ip(1:N/2);
        pam_symbols_temp2 = ip(N/2+1:N);
        pam_symbols_temp=pam_symbols_temp1+1i*pam_symbols_temp2;
        pam_symbols_te = rectpulse(pam_symbols_temp,Nsamp);
        Ase= transpose(lteZadoffChuSeq(1,nzc));
        Ase= rectpulse(Ase,Nsamp);
        pam_symbols_t = [pam_symbols_te(end-ncp+1:end) Ase    pam_symbols_te ];
        y0=[];
        fd=FD(ic);
        for i=1:length(pam_symbols_t)
            y0=[y0 (1*pam_symbols_t(i)*exp(1i*2*pi*(fc+fd)*[(i-1)*Tb+1/Fs:1/Fs:i*Tb]))];
        end
        
        
        dis=fas(ic);
        if(H==2 || H==3  )
            if(dis<1950)
                D1= D(1,ic);
            else
                D1= D(:,ic);
            end
        elseif(H==4)
            if(dis<1500)
                D1= D(1,ic);
            else
                D1= D(:,ic);
            end
        elseif(H==5)
            if(dis<1000)
                D1= D(1,ic);
            else
                D1= D(:,ic);
            end
        elseif(H==6)
            if(dis<2500)
                D1= D(1,ic);
            else
                D1= D(:,ic);
            end
            
        elseif(H==7)
            
            if(dis<2750)
                D1= D(1,ic);
            else
                D1= D(:,ic);
            end
            
        else
            D1=D(:,ic);
        end
        
        d0=D1;
        
        py=sum(y0*y0')/length(y0);
        y0=sqrt(1/py)*y0;
        y0=ncz(ic)*y0;
        
        
        for iu=1:length(d0)
            y(floor(d0(iu)*length(t)/Tb)+1:length(t)*(N/k*Nsamp+nzc*Nsamp+ncp)+floor((length(t)*d0(iu)/Tb)) )=...
                y(floor(d0(iu)*length(t)/Tb)+1:length(t)*(N/k*Nsamp+nzc*Nsamp+ncp)+floor((length(t)*d0(iu)/Tb)) )...
                -y0;
        end
    end
end
