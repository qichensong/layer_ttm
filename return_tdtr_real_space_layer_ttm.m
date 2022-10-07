function ft = return_tdtr_real_space_layer_ttm(td,omega_m,w0,w1,xx,yy,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
Ld =6/sqrt((w0^2+w1^2)/2);
nx1 = 80;
[qx,wq]=lgwt(nx1,0,Ld);
% qx = [qx(2)];
% wq = [1];
%wq = wq; % integrate from 0 to +infty

omegamax = 1e11*2*pi;
omega_rep = 1/(12.5e-9)*2*pi;
NN = floor(omegamax/omega_rep);
omega = (-NN:NN)*omega_rep+omega_m;
nw = length(omega);
nt = length(td);
irec(:,:,:) = freq_response_layer_ttm_iso(omega,w0,w1,qx,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2);
fr = zeros(size(xx,1),size(xx,2),nw,2);    

ii = 1;
jj = 1;
ff = numerical_int_hankel(irec,omega,qx,wq); % direct integration slower than
fr(ii,jj,:,1) =  abs(ff);
fp(:) = atan2(imag(ff),real(ff))*180/pi;
%          idx = find(fp>1e-5);
%          fp(idx) = fp(idx)-180;
fr(ii,jj,:,2) = fp;

ft = zeros(size(xx,1),size(xx,2),nt,2);
ft_temp = zeros(size(td));
for it = 1:length(td)
    ft_temp(it) = sum(ff.*exp(1j*omega*td(it)).*exp(-pi*(1e-10*omega/2/pi).^2));
end
ft(ii,jj,:,1) =  abs(ft_temp);
ft(ii,jj,:,2) = atan2(imag(ft_temp),real(ft_temp))*180/pi;

end