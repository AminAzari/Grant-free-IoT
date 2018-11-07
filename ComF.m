function co=ComF(MaxM)
coj=0;
MaxMI=MaxM>0;
SMaxM=sum(MaxMI,2);
t=size(MaxM,1);
if(t>10)
    error('error combi')
end
switch(t)
    case(1)
        for kk1=1:SMaxM(1)
            coj=coj+1;
            co(coj,1:t)=  [MaxM(1,kk1) ];
        end
    case(2)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                coj=coj+1;
                co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2)  ];
            end
        end
    case(3)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    coj=coj+1;
                    co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) ];
                end
            end
        end
    case(4)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        coj=coj+1;
                        co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4)];
                    end
                end
            end
        end
    case(5)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            coj=coj+1;
                            co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5)];
                        end
                    end
                end
            end
        end
        
    case(6)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            for kk6=1:SMaxM(6)
                                coj=coj+1;
                                co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5) MaxM(6,kk6)];
                            end
                        end
                    end
                end
            end
        end
    case(7)
        
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            for kk6=1:SMaxM(6)
                                for kk7=1:SMaxM(7)
                                    coj=coj+1;
                                    co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5) MaxM(6,kk6)  MaxM(7,kk7)];
                                end
                            end
                        end
                    end
                end
            end
        end
    case(8)
        
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            for kk6=1:SMaxM(6)
                                for kk7=1:SMaxM(7)
                                    for kk8=1:SMaxM(8)
                                        coj=coj+1;
                                        co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5) MaxM(6,kk6)  MaxM(7,kk7) MaxM(8,kk8)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    case(9)
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            for kk6=1:SMaxM(6)
                                for kk7=1:SMaxM(7)
                                    for kk8=1:SMaxM(8)
                                        for kk9=1:SMaxM(9)
                                            coj=coj+1;
                                            co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5) MaxM(6,kk6)  MaxM(7,kk7) MaxM(8,kk8) MaxM(9,kk9)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    case(10)
        
        for kk1=1:SMaxM(1)
            for kk2=1:SMaxM(2)
                for kk3=1:SMaxM(3)
                    for kk4=1:SMaxM(4)
                        for kk5=1:SMaxM(5)
                            for kk6=1:SMaxM(6)
                                for kk7=1:SMaxM(7)
                                    for kk8=1:SMaxM(8)
                                        for kk9=1:SMaxM(9)
                                            for kk10=1:SMaxM(10)
                                                coj=coj+1;
                                                co(coj,1:t)=  [MaxM(1,kk1) MaxM(2,kk2) MaxM(3,kk3) MaxM(4,kk4) MaxM(5,kk5) MaxM(6,kk6)  MaxM(7,kk7) MaxM(8,kk8) MaxM(9,kk9) MaxM(10,kk10)];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
end

