function [US]=SCALE(u,Xf)
% quando temos funções simétricas faz o respectivo preenchimento
Us=zeros(length(Xf));
    US=((Xf).*u.').';


    
    
