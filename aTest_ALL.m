%IHN
clc
clear all
close all
for Th=[3]

for La_Rg=[0.1:0.8:6]
    
%     
%     InSynT=1;
%     InSynF=1;
%     H=0;
%     for Nt=1:3
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     
%     InSynT=2;
%     InSynF=1;
%     H=0;
%     
%     for Nt=1:2
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     
%     
%     InSynT=2;
%     InSynF=2;
%     H=0;
%     
%     for Nt=1:2
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     
%     
%     
%     InSynT=1;
%     InSynF=1;
%     H=1;
%     for Nt=1:3
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     %
%     
%     InSynT=1;
%     InSynF=1;
%     H=2;
%     for Nt=[2,3,5]
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     
%     InSynT=1;
%     InSynF=1;
%     H=3;
%     for Nt=[2,3,5]
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     
%     
%     
%     InSynT=1;
%     InSynF=1;
%     H=4;
%     for Nt=[2,3,5]
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
%     %
    InSynT=1;
    InSynF=1;
    H=6;
    for Nt=[2,3,5]
        for lambda=La_Rg
            Test_final(lambda,Nt,InSynT,InSynF,H,Th);
        end
    end
%     
%     
%     InSynT=1;
%     InSynF=1;
%     H=6;
%     for Nt=[4]
%         for lambda=La_Rg
%             Test_final(lambda,Nt,InSynT,InSynF,H,Th);
%         end
%     end
    
    
    
    
end
end