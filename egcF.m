%IHN
function y=egcF(Au,cfo,y0,cy0,fc,Fs,Tb,ncp,nzc,Nsamp,Nt)
% t=[1/Fs:1/Fs:Tb];
% aak= FuncDecP_inf(t,y0,fc,Fs,Tb,cy0,ncp,nzc,Nsamp);
% rak=real(aak);
% Nl=length(aak);
% sy=sum(rak.^2)/Nl;
% 
% 
% for ij=1:Nt
%      yyk(ij,:)= FuncDecP_inf(t,Au(ij,:),fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
%      intf(ij,:)=real(aak)-real(yyk(ij,:));
%      sr(ij)=sum(intf(ij,:).^2)/Nl;
% end
% 
% sigr=sr;
% 
% B=diag(sigr)+sy*ones(Nt,Nt);
% C=sy*ones(Nt,1);
 A=ones(1,Nt);
y=zeros(size(Au(1,:)));
for i=1:Nt
   y=y+A(i)*Au(i,:);
end
y=y;%/(sqrt(sum(A)));

%  
%  yyk(1,:)= FuncDecP_inf(t,Au(1,:),fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
% yyk(2,:)= FuncDecP_inf(t,Au(2,:),fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
% 
% aak= FuncDecP_inf(t,y0,fc,Fs,Tb,cy0,ncp,nzc,Nsamp);
% 
% rak=real(aak);
% Nl=length(aak);
% sy=sum(rak.^2)/Nl;
% if1=real(aak)-real(yyk(1,:));
% if2=real(aak)-real(yyk(2,:));
% sr1=sum(if1.^2)/Nl;
% sr2=sum(if2.^2)/Nl;
% a1=sy/(sy*(1+sr1/sr2)+sr1);
% a2=sy/(sy*(1+sr2/sr1)+sr2);
% 
% sigr=[sr1,sr2];
% Nt=2;
% B=diag(sigr)+sy*ones(Nt,Nt);
% C=sy*ones(Nt,1);
% Ac=B\C;
% 
% y=(a1*Au(1,:)+a2*Au(2,:))/(a1+a2);
% 
%  
