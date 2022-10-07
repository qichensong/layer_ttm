function [te,tp,t2]=return_data_real_space_z_layer_ttm(omega,w0,w1,xx,yy,zz,zz2,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
Ld =6/sqrt((w0^2+w1^2)/2);
nx1 = 120;
[qx,wq]=lgwt(nx1,0,Ld);
% qx = [qx(2)];
% wq = [1];
%wq = wq; % integrate from 0 to +infty
nw = length(omega);
[tez,tpz,t2z] = freq_response_z_layer_ttm_iso(omega,zz,zz2,w0,w1,qx,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2);
fr = zeros(size(xx,1),size(xx,2),nw,2);    

ii = 1;
jj = 1;
te = zeros(length(zz),nw,2);
tp = zeros(length(zz),nw,2);
t2 = zeros(length(zz2),nw,2);
for iz = 1:length(zz)
    ff = numerical_int_hankel(reshape(tez(iz,:,:),[nw,nx1]),omega,qx,wq); % direct integration slower than
    te(iz,:,1) =  abs(ff);
    te(iz,:,2) = atan2(imag(ff),real(ff))*180/pi;
    ff = numerical_int_hankel(reshape(tpz(iz,:,:),[nw,nx1]),omega,qx,wq); % direct integration slower than
    tp(iz,:,1) =  abs(ff);
    tp(iz,:,2) = atan2(imag(ff),real(ff))*180/pi;
end
for iz = 1:length(zz2)
    ff = numerical_int_hankel(reshape(t2z(iz,:,:),[nw,nx1]),omega,qx,wq); % direct integration slower than
    t2(iz,:,1) =  abs(ff);
    t2(iz,:,2) = atan2(imag(ff),real(ff))*180/pi;
end
end