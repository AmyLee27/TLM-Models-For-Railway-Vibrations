function [feMat] = trackMat(feElt,Ezz_s,nu_s,rho_s,D_s,Ezz_a,nu_a,rho_a,D_a,Ezz_b,nu_b,rho_b,D_b,Ezz_usl,nu_usl,rho_usl,D_usl,Ezz_sl,nu_sl,rho_sl,D_sl,Ezz_rp,nu_rp,rho_rp,D_rp,E_r,nu_r,rho_r,D_r)

feMat = zeros(numel(feElt(:,1)),5);

feMat(:,1:2) = feElt(:,1:2);

index1 = feElt(feElt(:,2)==1,1);
index2 = feElt(feElt(:,2)==2,1);
index3 = feElt(feElt(:,2)==3,1);
index4 = feElt(feElt(:,2)==4,1);
index5 = feElt(feElt(:,2)==5,1);
index6 = feElt(feElt(:,2)==6,1);
index7 = feElt(feElt(:,2)==7,1);

feMat(index1,3:6) = repmat([Ezz_s,nu_s,rho_s,D_s],numel(index1),1);
feMat(index2,3:6) = repmat([Ezz_a,nu_a,rho_a,D_a],numel(index2),1);
feMat(index3,3:6) = repmat([Ezz_b,nu_b,rho_b,D_b],numel(index3),1);
feMat(index4,3:6) = repmat([Ezz_usl,nu_usl,rho_usl,D_usl],numel(index4),1);
feMat(index5,3:6) = repmat([Ezz_sl,nu_sl,rho_sl,D_sl],numel(index5),1);
feMat(index6,3:6) = repmat([Ezz_rp,nu_rp,rho_rp,D_rp],numel(index6),1);
feMat(index7,3:6) = repmat([E_r,nu_r,rho_r,D_r],numel(index7),1);

end