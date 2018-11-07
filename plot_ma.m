
%IHN
%IHN
clc
clear all
close all
Pcl=0.080;
Pt=10^((21-30)/10);
alf=1.1 ;
Pc=0.055;
Leng_pack= 0.5;
La_Rg=.1:.4:4;
% lambda=La_Rg;
Fm=200;
N=50;%bits/pack
% Leng_pack=0.5;%%%%%%%%%%%%   Eth  Mohem
Es=5*Pc;
E0=5000;
Tw=5;
Rc1=50;
Rc2=3000;
Esyt=Pcl*3;
Dsyt=2;
Esytf=Pcl*3;
Dsytf=3;
FFh= [100:250:3000];
nFFh= [100:100:3000];
Nd=20000;
FasD = sqrt((Rc2^2-Rc1^2)*rand(Nd,1)+Rc1^2);
Ind=zeros(Nd,length(nFFh));
for io=1:length(nFFh)
    Ind(:,io)=  (nFFh(io)-50<FasD & FasD<nFFh(io)+50);
    sIn(io)=sum( Ind(:,io));
end
% save('dataI.mat','sIn')
% load dataI.mat
Th=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
InSynT=1;
InSynF=1;
Esy=0;Dsy=0;
H=0;
for Th=[1]
    
    for kl=1
        if(Th==1 && kl==2)
            
            InSynT=2;
            InSynF=1;
            Esy=Esyt;Dsy=Dsyt;
        end
        if(Th==1 && kl==3)
            InSynT=2;
            InSynF=2;
            Esy=Esytf;Dsy=Dsytf;
        end
        
        if(Th==1 && kl==2)
            
            InSynT=2;
            InSynF=1;
            Esy=Esyt;Dsy=Dsyt;
        end
        if(Th==1 && kl==1)
            InSynT=1;
            InSynF=1;
            Esy=0;Dsy=0;
        end
        
        
        Nt=1;
        co=0;
        for lambda=La_Rg
            co=co+1;
            X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th)  ];
            load([X,'.mat'])
            %     SR111(co)=mean(SucRate);
            NSR111(co,:)=nanmean(nNSucRate);
        end
        
        EE=0;
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR).*Pc*Tw+Esy) ;
        end
        
        figure(1)
        plot(La_Rg ,EE/1000)
        hold on
        
        DD=0;
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
        end
        figure(2)
        plot(La_Rg ,DD)
        hold on
        
        LT=0;sLT=zeros(length(nFFh),length(La_Rg));
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR).*Pc*Tw+Esy);
        end
        LiF=[];
        for ct=1:length(La_Rg)
            cot=1;
            for cl=1:length(nFFh)
                LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
                cot=cot+sIn(cl);
            end
        end
        for ok=1:length(La_Rg)
            JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
        end
        % figure(9)
        % plot(La_Rg ,LT2)
        % hold on
        figure(3)
        plot(La_Rg ,LT)
        hold on
        
        
        figure(8)
        plot(La_Rg ,JF)
        hold on
        
        
        kSR=0;
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            kSR=kSR+(sIn(j)/Nd)*SR;
        end
        
        figure(4)
        plot(La_Rg ,kSR)
        hold on
        
        
        kSR1=0;
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            kSR1=kSR1+(sIn(j))*SR;
            sSR(1,j)=SR(end);
            sSR(2,j)=sIn(j)/Nd;
        end
        fkSR1=kSR1.^2/Nd;
        
        kSR2=0;
        for j=1:length(nFFh)
            SR=NSR111(:,j)';
            kSR2=kSR2+(sIn(j))*SR.^2;
        end
        fkSR2=kSR2;
        
        figure(5)
        plot(La_Rg ,fkSR1./fkSR2)
        hold on
        
        co=0;
        lambda=2.1;
        for n=FFh
            co=co+1;
            X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
            load([X,'.mat'])
            ANSR(co)=nanmean(NSucRate(:,co));
        end
        figure(6)
        plot(FFh ,ANSR)
        hold on
        
        
        %%%--------------------------
        
        Nt=2;
        co=0;
        for lambda=La_Rg
            co=co+1;
            X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
            load([X,'.mat'])
            %     SR111(co)=mean(SucRate);
            NSR211(co,:)=nanmean(nNSucRate);
        end
        
        
        EE=0;
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(2*Leng_pack*(alf*Pt +Pc)+Pc*Leng_pack)+(1./SR).*Pc*Tw+Esy) ;
        end
        
        figure(1)
        plot(La_Rg ,EE/1000)
        hold on
        
        
        DD=0;
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            NSR=NSR111(:,j)';
            Ps=(1-(1-NSR).^2);
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(3*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack);
        end
        figure(2)
        plot(La_Rg ,DD)
        hold on
        
        
        LT=0;sLT=zeros(length(nFFh),length(La_Rg));
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*(2*Leng_pack*(alf*Pt+Pc )+Pc*Leng_pack)+(1./SR).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*(2*Leng_pack*(alf*Pt+Pc )+Pc*Leng_pack)+(1./SR).*Pc*Tw+Esy);
        end
        
        figure(3)
        plot(La_Rg ,LT)
        hold on
        
        
        LiF=[];
        for ct=1:length(La_Rg)
            cot=1;
            for cl=1:length(nFFh)
                LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
                cot=cot+sIn(cl);
            end
        end
        for ok=1:length(La_Rg)
            JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
        end
        figure(8)
        plot(La_Rg ,JF)
        hold on
        
        kSR=0;
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            kSR=kSR+(sIn(j)/Nd)*SR;
        end
        
        figure(4)
        plot(La_Rg ,kSR)
        hold on
        
        kSR1=0;
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            kSR1=kSR1+(sIn(j))*SR;
        end
        fkSR1=kSR1.^2/Nd;
        
        kSR2=0;
        for j=1:length(nFFh)
            SR=NSR211(:,j)';
            kSR2=kSR2+(sIn(j))*SR.^2;
        end
        fkSR2=kSR2;
        
        figure(5)
        plot(La_Rg ,fkSR1./fkSR2)
        hold on
        
        
        
        co=0;
        lambda=2.1;
        for n=FFh
            co=co+1;
            X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
            load([X,'.mat'])
            ANSR(co)=nanmean(NSucRate(:,co));
        end
        figure(6)
        plot(FFh ,ANSR)
        hold on
        %%%--------------------------
    end
    
    InSynT=1;
    InSynF=1;
    Esy=0;Dsy=0;
    
    
    Nt=3;
    co=0;
    for lambda=La_Rg
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        %     SR111(co)=mean(SucRate);
        NSR311(co,:)=nanmean(nNSucRate);
    end
    
    
    EE=0;
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(3*Leng_pack*(alf*Pt+Pc )+2*Pc*Leng_pack)+(1./SR).*Pc*Tw+Esy) ;
    end
    
    figure(1)
    plot(La_Rg ,EE/1000)
    hold on
    
    DD=0;
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        NSR=NSR111(:,j)';
        Ps=(1-(1-NSR).^3);
        DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
    end
    figure(2)
    plot(La_Rg ,DD)
    hold on
    
    
    
    LT=0;sLT=zeros(length(nFFh),length(La_Rg));
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt +Pc)+2*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
        sLT(j,:)=E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt +Pc)+2*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
    end
    
    figure(3)
    plot(La_Rg ,LT)
    hold on
    LiF=[];
    for ct=1:length(La_Rg)
        cot=1;
        for cl=1:length(nFFh)
            LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
            cot=cot+sIn(cl);
        end
    end
    for ok=1:length(La_Rg)
        JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
    end
    figure(8)
    plot(La_Rg ,JF)
    hold on
    
    kSR=0;
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        kSR=kSR+(sIn(j)/Nd)*SR;
    end
    
    figure(4)
    plot(La_Rg ,kSR)
    hold on
    
    
    kSR1=0;
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        kSR1=kSR1+(sIn(j))*SR;
    end
    fkSR1=kSR1.^2/Nd;
    
    kSR2=0;
    for j=1:length(nFFh)
        SR=NSR311(:,j)';
        kSR2=kSR2+(sIn(j))*SR.^2;
    end
    fkSR2=kSR2;
    
    figure(5)
    plot(La_Rg ,fkSR1./fkSR2)
    hold on
    
    co=0;
    lambda=2.1;
    for n=FFh
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        ANSR(co)=nanmean(NSucRate(:,co));
    end
    figure(6)
    plot(FFh ,ANSR)
    hold on
end
%%%--------------------------

%%%--------------------------
for gf=1:1
Nres=10;
Tpr=2;
Nt=1;
co=0;
Esy=Esytf;Dsy=Dsytf;
for lambda=La_Rg
    co=co+1;
    %     load datap.mat
    Pss(co)=FunPs(Nres,lambda,Tpr);
    if(1/(11/Tpr-lambda) <0)
        Del(co)=nan;
    else
        Del(co)=Tpr/2+min(2*Tpr,1/(11/Tpr-lambda));
    end
end

lambda=0.1;
X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
load([X,'.mat'])
SSRn=nanmean(nNSucRate);



figure(1)
EE=0;
for j=1:length(nFFh)
    SRn=SSRn(j);
    EE=EE+(sIn(j)/Nd)*La_Rg.*N./(Es+Esy+Pc*Del+(1./SRn).*(1./Pss).*(Leng_pack/3*(alf*Pt+Pc)+Pcl*1)+(1./SRn)*(Leng_pack*(alf*Pt+Pc)+Pcl*1+Pc*Tw));
end
plot(La_Rg ,EE/1000)
hold on

DD=0;
for j=1:length(nFFh)
    SRn=SSRn(j);
    DD=DD+(sIn(j)/Nd)*min(20,(Dsy+Del +(1./Pss-1).*Tpr+(1./SRn)*(Leng_pack)+(1./SRn-1)*Tw ));
end
figure(2)
plot(La_Rg ,DD)
hold on

LT=0;sLT=zeros(length(nFFh),length(La_Rg));

for j=1:length(nFFh)
    SRn=SSRn(j);
    LT=LT+(sIn(j)/Nd)*(E0./(Es+Esy+Pc*Del+(1./SRn).*(1./Pss).*(Leng_pack/3*(alf*Pt+Pc)+Pcl*1)+(1./SRn)*(Leng_pack*(alf*Pt+Pc)+Pcl*1+Pc*Tw)));
    sLT(j,:)=(E0./(Es+Esy+Pc*Del+(1./SRn).*(1./Pss).*(Leng_pack/3*(alf*Pt+Pc)+Pcl*1)+(1./SRn)*(Leng_pack*(alf*Pt+Pc)+Pcl*1+Pc*Tw)));
end

figure(3)
plot(La_Rg ,LT)
hold on

LiF=[];
for ct=1:length(La_Rg)
    cot=1;
    for cl=1:length(nFFh)
        LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
        cot=cot+sIn(cl);
    end
end
for ok=1:length(La_Rg)
    JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
end
figure(8)
plot(La_Rg ,JF)
hold on 

% figure(9)
% plot(La_Rg ,sLT)
% hold on

kSR=0;
for j=1:length(nFFh)
    SRn=SSRn(j);
    kSR=kSR+(sIn(j)/Nd)*SRn*Pss;
end

figure(4)
plot(La_Rg ,kSR)
hold on


kSR1=0;
for j=1:length(nFFh)
    SRn=SSRn(j);
    kSR1=kSR1+(sIn(j))*SRn*Pss;
end
fkSR1=kSR1.^2/Nd;

kSR2=0;
for j=1:length(nFFh)
    SRn=SSRn(j);
    kSR2=kSR2+(sIn(j))*(SRn*Pss).^2;
end
fkSR2=kSR2;

figure(5)
plot(La_Rg ,fkSR1./fkSR2)
hold on

%         co=0;
% lambda=2.1;
%         for n=FFh
%             co=co+1;
%             X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
%             load([X,'.mat'])
%             ANSR(co)=nanmean(NSucRate(:,co));
%         end

lambda=0.1;
X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
load([X,'.mat'])
nSRn=nanmean(NSucRate);
nPss=FunPs(Nres,lambda,Tpr);
figure(6)
plot(FFh ,Pss(6)*nSRn)
hold on

lambda=5.7;
X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
load([X,'.mat'])
sSR(3,:)=nanmean(nNSucRate);
end

%%%--------------------------

% % %
% % % % % % % %%%--------------------------
Th=1;
Esy=0;Dsy=0;
for H=[2 ,4,5]
    switch(H)
        case(5)
            UB=11;
        case(4)
            UB=15;
        case(2)
            UB=19;
    end
%     H=2;
    Nt=3;
    co=0;
    for lambda=La_Rg
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        SR111(co)=mean(SucRate);
        NSR4113(co,:)=nanmean(nNSucRate);
    end
    
    
    EE=0;
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        if(j<UB)
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR ).*Pc*Tw+Esy) ;
        else
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(3*Leng_pack*(alf*Pt+Pc )+2*Pc*Leng_pack)+(1./SR).*Pc*Tw+Esy) ;
        end
    end
    
    figure(1)
    plot(La_Rg ,EE/1000)
    hold on
    
    DD=0;
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        NSR=NSR111(:,j)';
        if(j<UB)
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
            
        else
            Ps=(1-(1-NSR).^3);
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
        end
    end
    figure(2)
    plot(La_Rg ,DD)
    hold on
    
    
    
    LT=0;sLT=zeros(length(nFFh),length(La_Rg));
    
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        if(j<UB)
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc )+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc )+(1./SR ).*Pc*Tw+Esy);
        else
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt+Pc )+2*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt+Pc )+2*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
        end
    end
    
    figure(3)
    plot(La_Rg ,LT)
    hold on
    % %
    LiF=[];
    for ct=1:length(La_Rg)
        cot=1;
        for cl=1:length(nFFh)
            LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
            cot=cot+sIn(cl);
        end
    end
    for ok=1:length(La_Rg)
        JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
    end
    figure(8)
    plot(La_Rg ,JF)
    hold on
    
    kSR=0;
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        kSR=kSR+(sIn(j)/Nd)*SR;
    end
    
    figure(4)
    plot(La_Rg ,kSR)
    hold on
    
    kSR1=0;
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        kSR1=kSR1+(sIn(j))*SR;
    end
    fkSR1=kSR1.^2/Nd;
    
    kSR2=0;
    for j=1:length(nFFh)
        SR=NSR4113(:,j)';
        kSR2=kSR2+(sIn(j))*SR.^2;
    end
    fkSR2=kSR2;
    
    figure(5)
    plot(La_Rg ,fkSR1./fkSR2)
    hold on
    
    co=0;
    lambda=2.1;
    for n=FFh
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        ANSR(co)=nanmean(NSucRate(:,co));
    end
    figure(6)
    plot(FFh ,ANSR)
    hold on
    %%--------------------------
    
    %%--------------------------
    %%--------------------------
%     H=2;
    Nt=5;
    co=0;
    for lambda=La_Rg
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        SR111(co)=mean(SucRate);
        NSR4115(co,:)=nanmean(nNSucRate);
    end
    
    
    EE=0;
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        if(j<UB)
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR ).*Pc*Tw+Esy) ;
        else
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(5*Leng_pack*(alf*Pt  +Pc)+4*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy) ;
        end
    end
    
    figure(1)
    plot(La_Rg ,EE/1000)
    hold on
    
    DD=0;
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        NSR=NSR111(:,j)';
        if(j<UB)
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
            
        else
            Ps=(1-(1-NSR).^3);
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
        end
    end
    figure(2)
    plot(La_Rg ,DD)
    hold on
    
    
    
    LT=0;sLT=zeros(length(nFFh),length(La_Rg));
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        if(j<UB)
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR ).*Pc*Tw+Esy);
        else
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*(5*Leng_pack*(alf*Pt +Pc )+4*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*(5*Leng_pack*(alf*Pt +Pc )+4*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
        end
    end
    
    figure(3)
    plot(La_Rg ,LT)
    hold on
    
    % %
    LiF=[];
    for ct=1:length(La_Rg)
        cot=1;
        for cl=1:length(nFFh)
            LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
            cot=cot+sIn(cl);
        end
    end
    for ok=1:length(La_Rg)
        JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
    end
    figure(8)
    plot(La_Rg ,JF)
    hold on% % %
    
    
    kSR=0;
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        kSR=kSR+(sIn(j)/Nd)*SR;
    end
    
    figure(4)
    plot(La_Rg ,kSR)
    hold on
    
    kSR1=0;
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        kSR1=kSR1+(sIn(j))*SR;
    end
    fkSR1=kSR1.^2/Nd;
    
    kSR2=0;
    for j=1:length(nFFh)
        SR=NSR4115(:,j)';
        kSR2=kSR2+(sIn(j))*SR.^2;
    end
    fkSR2=kSR2;
    
    figure(5)
    plot(La_Rg ,fkSR1./fkSR2)
    hold on
    
    co=0;
    lambda=2.1;
    for n=FFh
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        ANSR(co)=nanmean(NSucRate(:,co));
    end
    figure(6)
    plot(FFh ,ANSR)
    hold on
    %%--------------------------
    
    %%--------------------------
%     H=2;
    Nt=7;
    co=0;
    for lambda=La_Rg
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        SR111(co)=mean(SucRate);
        NSR4117(co,:)=nanmean(nNSucRate);
    end
    
    
    EE=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        if(j<UB)
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR ).*Pc*Tw+Esy) ;
        else
            EE=EE+(sIn(j)/Nd)*La_Rg.*N./((1./SR).*(7*Leng_pack*(alf*Pt  +Pc)+6*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy) ;
        end
    end
    
    figure(1)
    plot(La_Rg ,EE/1000)
    hold on
    
    DD=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        NSR=NSR111(:,j)';
        if(j<UB)
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
            
        else
            Ps=(1-(1-NSR).^7);
            DD=DD+(sIn(j)/Nd)*min(20,Dsy+(1./SR-1).*(9*Leng_pack+Tw)+...
                NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack+...
                +(1-NSR).^3.*(NSR)./Ps.*7*Leng_pack+(1-NSR).^4.*(NSR)./Ps.*9*Leng_pack);
        end
    end
    figure(2)
    plot(La_Rg ,DD)
    hold on
    
    
    
    LT=0;sLT=zeros(length(nFFh),length(La_Rg));
    LT2=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        if(j<UB)
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR ).*Pc*Tw+Esy);
        else
            LT=LT+(sIn(j)/Nd)*E0./(Es+(1./SR).*(7*Leng_pack*(alf*Pt +Pc )+6*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
            sLT(j,:)=E0./(Es+(1./SR).*(7*Leng_pack*(alf*Pt +Pc )+6*Pc*Leng_pack)+(1./SR ).*Pc*Tw+Esy);
        end
    end
    
    figure(3)
    plot(La_Rg ,LT)
    hold on
    
    % %
    LiF=[];
    for ct=1:length(La_Rg)
        cot=1;
        for cl=1:length(nFFh)
            LiF(ct,cot:cot+sIn(cl)-1)=sLT(cl,ct);
            cot=cot+sIn(cl);
        end
    end
    for ok=1:length(La_Rg)
        JF(ok)=(sum(LiF(ok,:))^2/20000)/sum(LiF(ok,:).^2);
    end
    figure(8)
    plot(La_Rg ,JF)
    hold on% % %
    kSR=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        kSR=kSR+(sIn(j)/Nd)*SR;
    end
    
    figure(4)
    plot(La_Rg ,kSR)
    hold on
    
    kSR1=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        kSR1=kSR1+(sIn(j))*SR;
    end
    fkSR1=kSR1.^2/Nd;
    
    kSR2=0;
    for j=1:length(nFFh)
        SR=NSR4117(:,j)';
        kSR2=kSR2+(sIn(j))*SR.^2;
    end
    fkSR2=kSR2;
    
    figure(5)
    plot(La_Rg ,fkSR1./fkSR2)
    hold on
    
    co=0;
    lambda=2.1;
    for n=FFh
        co=co+1;
        X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
        load([X,'.mat'])
        ANSR(co)=nanmean(NSucRate(:,co));
    end
    figure(6)
    plot(FFh ,ANSR)
    hold on
end

%%---------------
%,'Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7'
%    '1950-Hy3','1950-Hy5','1950-Hy7','1500-Hy3','1500-Hy5','1500-Hy7','1000-Hy3','1000-Hy5','1000-Hy7'

figure(1)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('Energy Efficiency (Bit/mJoule)')

figure(2)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('Service Delay (Sec)')

figure(3)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('Average Battery Lifetime (\times reporting period)')

figure(4)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('Probability of success in transmission')

figure(5)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('fairness  ')

figure(6)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('Ac-vs-dis  ')

figure(8)
legend('N=1; TiSy; FrSy','N=2; TiSy; FrAs','N=3; TiAs; FrAs','Grant-Based','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7','Hy3','Hy5','Hy7')
grid on
xlabel('\Lambda: Aggregated arrival rate of packets')
ylabel('fairness-lif  ')
% % %% ------------
% % % % %%
% % close all
% % H=0;
% %
% % InSynT=1;
% % InSynF=1;
% % fi=0;
% % for lambda=La_Rg(1)
% %     fi=fi+1;
% %     for Nt=[2:3]
% %
% %         co=0;
% %         %     lambda=2.1;
% %         for n=FFh
% %             co=co+1;
% %             X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H) ,'Th',num2str(Th) ];
% %             load([X,'.mat'])
% %             ANSR(co)=nanmean(NSucRate(:,co));
% %         end
% %         SR=  ANSR;
% %         ko=0;
% %         for dis =FFh
% %             ko=ko+1;
% %             % dis=nFFh(5);
% %             Plu=133 + 38.3*log(dis/1000);
% %             W=200;
% %
% %             Pot=21-30;
% %             PiGz=(10.^((Pot-Plu)/10));
% %             No=10^(-20.4)*W;
% %             del=3.83;
% %             Ra=3000;
% %
% %
% %
% %             Th=1;
% %             Pn =exp(-No*Th/PiGz );
% %             lam= Nt*(W/(2*Fm+W))*(1/(pi*Ra^2))*lambda*Leng_pack;
% %             Ca=(2*(pi^2))/del*csc(2*pi/del); %pi/sinc(2/del)
% %             Pp=exp(-lam*(dis^2)*((Th)^(2/del))*Ca);
% %             PS(ko)=Pn*Pp;
% %             Pnn(ko)=Pn;
% %             Ppp(ko)=Pp;
% %             fun=@(x)exp(-No*Th./(x*PiGz)).*exp(-lam.*(dis.^2).*((Th./x).^(2/del)).*Ca).*x.*exp(-x);
% %             Psmr(ko)=integral(fun,0,10000);
% %
% % %             Thm=Th/Nt;
% % %             Pnm =exp(-No*Thm/PiGz );
% % %             lam= Nt*(1/(pi*Ra^2))*(W/(2*Fm+W))*lambda*Leng_pack;
% % %             Ca=(2*(pi^2))/del*csc(2*pi/del); %pi/sinc(2/del)
% % %             Ppm=exp(-lam*(dis^2)*(Thm^(2/del))*Ca);
% % %             PSm(ko)=Pnm*Ppm;
% %         end
% %         figure(fi)
% %
% %
% %         plot(FFh,PS,'-b')
% %         hold on
% %         plot(FFh,Ppp,'-c')
% %
% %         % plot(La_Rg,Pp,'-r')
% %         plot(FFh,1-(1-PS).^Nt,'-g')
% %
% % %         plot(FFh,Psmr,'-r')
% %         plot(FFh,SR,'-m')
% %          plot(FFh,Pnn,'-k')
% %         legend('PS1','ppp','PSNt','PSSim','Pn')
% % %     Ppp
% %     end
% % end
% % % hold off
% %
