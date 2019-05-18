function [Kdy]=Kdyeleinf1(prop,k1,k2,w)
% Calcula a matriz Kdyn do elemento
% Calculates the Kdyn matrix element

% - x - matriz [zi], com as coordenadas dos nós
%                    with the coordinates of the nodes
% - ro - massa volúmica
%        density
% - k - número de onda (rad/m)
%       wave number
% - E - módulo de Young;
%       Young's modulus
% - niu - coeficiente de Poisson;
%         Poisson's ratio
% - del - amortecimento;
%         damping
% - w - frequencia (rad/s)
%       frequency
% - prop- vector (linha, com as propriedades do elemento)
%         vector ( line with the element properties )
Eh=prop(1);
niuh=prop(2);
delh=prop(3);
mh=prop(4);

          if k1==0,
             fi=pi/2;
          else
             fi=abs(atan(k2/k1));
          end

k3=(k1^2+k2^2)^0.5;

la=niuh*Eh*(1+2*i*delh*sign(w))/((1+niuh)*(1-2*niuh));% attention - confirm the signal - formulation of takemyia
mi=Eh*(1+2*i*delh*sign(w))/(2*(1+niuh));% attention - confirm the signal - formulation of takemiya
cl=((la+2*mi)/mh)^0.5;
ct=(mi/mh)^0.5;
al=(k3^2-w^2/cl^2)^0.5; %attention, must be confirmed 
at=(k3^2-w^2/ct^2)^0.5; %attention, must be confirmed
bl=w^2/cl^2;

if w==0
    R=[0 1 0;0 0 1;-(la+3*mi)/(2*al*mi) 0 i*k3/at];
    S=[0 -mi*at 0;-i*k3*mi/al 0 -1*mi*(k3^2+at^2)/at;la+2*mi 0 -2*i*mi*k3];
else
    
R=[0 1 0; -i*k3/bl 0 1; al/bl 0 i*k3/at];
S=[0 -mi*at 0; 2*i*mi*al*k3/bl 0 -mi*(k3^2+at^2)/at;la-2*mi*al^2/bl 0 -2*i*mi*k3];
   
end


Qpo=-R*inv(S);
Qpi=[sin(fi)*sign(k2),cos(fi)*sign(k1),0;-cos(fi)*sign(k1),sin(fi)*sign(k2),0;0,0,1]*Qpo*[sin(fi)*sign(k2),-cos(fi)*sign(k1),0;cos(fi)*sign(k1),sin(fi)*sign(k2),0;0,0,1];


Kdy=inv(Qpi);
