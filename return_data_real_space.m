function fr = return_data_real_space(omega,w0,w1,xx,yy,rhocp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
Ld =6/sqrt((w0^2+w1^2)/2);
nx1 = 120;
[qx,wq]=lgwt(nx1,0,Ld);
% qx = [qx(2)];
% wq = [1];
%wq = wq; % integrate from 0 to +infty
nw = length(omega);
rho_cp = [rhocp+Ce,0,C2];
Sr = [kzz+kez,Gpp,k2z];
Sz = [kzz+kez,Gpp,k2z];
d = [L(1),0,L(2)];
alpha = alp;
sep = 0;

temp = freq_response(omega,w0,w1,sep,rho_cp,Sr,Sz,d,alpha);
fr(:,1) = abs(temp);
fr(:,2) = atan2(imag(temp),real(temp))*180/pi;

end