function [stre1,stre2,stre3,stre4,stre5,stre6]=Stress_p(strt1,strt2,strt3,strt4,strt5,strt6,properties1,w,elestrain);
for j=1:length(elestrain);
    [D]=Dmatrix(properties1(elestrain(j),1),properties1(elestrain(j),2),properties1(elestrain(j),3),w);
    strai=[strt1(j);strt2(j);strt3(j);strt4(j);strt5(j);strt6(j)];
    stree=D*strai;
    stre1(j)=stree(1,1);
    stre2(j)=stree(2,1);
    stre3(j)=stree(3,1);
    stre4(j)=stree(4,1);
    stre5(j)=stree(5,1);
    stre6(j)=stree(6,1);

end

