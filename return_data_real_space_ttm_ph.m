function fr = return_data_real_space_ttm_ph(omega,w0,w1,xx,yy,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g)
Ld = 6/sqrt((w0^2+w1^2)/2);
nx1 = 80;
[qx,wq]=lgwt(nx1,-Ld,Ld);
%wq = wq; % integrate from -infty to +infty
nw = length(omega);
irec(:,:,:) = freq_response_tensor_qgrid_ttm_ph(omega,w0,w1,qx,rho_cp,kzz,kxx,alp,gamma,Ce,kez,kex,g);
fr = zeros(size(xx,1),size(xx,2),nw,2);    

for ii = 1:size(xx,1)
    for jj = 1:size(xx,2)
         ff = numerical_int(irec,omega,qx,wq,xx(ii,jj),yy(ii,jj)); % direct integration slower than
         fr(ii,jj,:,1) =  abs(ff);
         fp(:) = atan2(imag(ff),real(ff))*180/pi;
%          idx = find(fp>1e-5);
%          fp(idx) = fp(idx)-360;
         fr(ii,jj,:,2) = fp;
    end
end
end