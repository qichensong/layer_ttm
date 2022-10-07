function fr = return_data_real_space_layer_ttm(omega,w0,w1,xx,yy,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
Ld =6/sqrt((w0^2+w1^2)/2);
nx1 = 120;
[qx,wq]=lgwt(nx1,0,Ld);
% qx = [qx(2)];
% wq = [1];
%wq = wq; % integrate from 0 to +infty
nw = length(omega);
irec(:,:) = freq_response_layer_ttm_iso(omega,w0,w1,qx,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2);
fr = zeros(size(xx,1),size(xx,2),nw,2);    

ii = 1;
jj = 1;
ff = numerical_int_hankel(irec,omega,qx,wq); % direct integration slower than
fr(ii,jj,:,1) =  abs(ff);
fp(:) = atan2(imag(ff),real(ff))*180/pi;
%          idx = find(fp>1e-5);
%          fp(idx) = fp(idx)-180;
fr(ii,jj,:,2) = fp;

end