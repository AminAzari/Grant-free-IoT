%% -TRansmitter%%%%Baseband%%%%%
function [yo,ip,y0,cy0,fas,zar,ncz]=FunTxP(N,k,Nsamp, ncp,nzc,fc,Fs,Tb,FD,D,cons,Pt,Rc1,Rc2,W,Nt,H)
t=[1/Fs:1/Fs:Tb];
dmx=max(D);
dmax=max(dmx);
Ntx1=dmax+2*(N/k*Nsamp+ncp+nzc*Nsamp)*Tb;
Ntx=floor(Ntx1* length(t))/( length(t));

ip=zeros(size(D,1),N);
yo=0;
zar=zeros(Nt,length(FD));
for i=1:length(FD)
    dis = sqrt((Rc2^2-Rc1^2)*rand(1,1)+Rc1^2);
    if(H==0)
        pH=1; 
        D1=D(:,i);
    elseif(H==1)
        pH=Nt;
        D1=D(:,i);
    elseif(H==2)
        pH=1;
        if(dis<1950)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
    elseif(H==3)
        if(dis<1950)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
        pH=length(D1);
    elseif(H==4)
        pH=1;
        if(dis<1500)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
    elseif(H==5)
        pH=1;
        if(dis<1000)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
         elseif(H==6)
        pH=1;
        if(dis<2500)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
        
        elseif(H==7)
        pH=1;
        if(dis<2750)
            D1= D(1,i);
        else
            D1= D(:,i);
        end
        
    end
    
    [yy1,ip(i,:),y0,cy0,zar1]=FunGen(N, Nsamp,Fs,Tb,fc,FD(i),D1,ncp,nzc,Ntx,k,cons,dis,Pt,pH);
    fas(i)=dis;
    PLr=133 + 38.3*log(dis./1000);
    %     zar(i)=sqrt((10^(Pt/(10)))/(10^(PLr/10))/pH );
    %     yy2=yy1*zar(i);
    ncz(i,1)=sqrt((10^(Pt/(10)))/(10^(PLr/10)) );
    %     n=1/sqrt(2)*[randn(1,length(yy1)) + 1i*randn(1,length(yy1))];
    %     No=sqrt(10^(-20.4)*W);
    
    yy3 = yy1 ;%+No*n;  %additive white gaussian noise%+awgn(y,snr_rg(ii));
    yo=yo+  yy3;
    
    zar(1:length(D1),i)=zar1;
    
end



% % 
% % 
% % periodogram(yo);
% % AxHandle=gca;
% % H = findobj(AxHandle, '-property', 'YData');
% % HX=H.XData;
% % HY=H.YData;
% % [~, po]=findpeaks(HY,'MinPeakDistance',10,'MinPeakHeight',25);
% % round(HX(po)*Fs/2)
% % 2
% % 
% % 
% % 
% % 
% % 
% % %% freq analysis of tx signal
% % %
% % figure()
% % y_rj=yo;
% % FX=[1:1:Fs-1];
% % periodogram(y_rj(1:Fs*Tb),[],FX,Fs);
% % figure()
% % periodogram(y_rj(1:Nsamp*Fs*Tb),[],Fs,'centered');
% % % AxHandle=gca;
% % % for H = findobj(AxHandle, '-property', 'YData')
% % %   set(H, 'YData', 10.^((get(H,'YData'))/10));
% % % end
% % 
% % figure()
% % periodogram(y_rj(1:Nsamp*Fs*Tb));
% % figure()
% % periodogram(y_rj(1:Nsamp*Fs*Tb),[],FX,Fs);
% % 
% % figure()
% % periodogram(y_rj,[],FX,Fs);
% % 2
% % 
% % 
% % % az koja dige resolution kafi hast: az 7
% % maxd=10;
% % for d=7:maxd
% %     X=yo(1:d*Nsamp*Fs*Tb);
% %     figure()
% %     periodogram(X,[],Fs,'centered');
% % end
% % 2

