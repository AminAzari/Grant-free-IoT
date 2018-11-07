%IHN
function rkk=FunRaJ(x ,M,ord)

rkk=1;
ml=length(M);
for i=1:ord 
    p1=x+i;
    p2=x-i;
    if(p1>ml)
        p1=ml;
    end
    if(p2<1)
        p2=1; 
    end
    rkk=rkk*(M(x)/M(p1)*M(x)/M(p2));
end