function y = dif_te(x,param_fit,omega,amp,phase,param0,power,fit_type)
param = param0;
for i = 1:length(param_fit)
    param(param_fit(i)) = x(i);
end
% param = [w0,C,kz,kz,alp,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR];
w0 = param(1);
C = param(2);
kx = param(3);
kz = param(4);
alp = param(5);
Ce = param(6);
kex = param(7);
kez = param(8);
g = param(9);
L = [param(10),param(11)];
Gep = param(12);
Gpp = param(13);
k2x = param(14);
k2z = param(15);
C2 = param(16);
CTR = param(17);

data = return_data_real_space_layer_ttm(omega,w0,w0,[0],[0],C,kz,kz,alp,alp,Ce,kez,kez,g,L,Gep,Gpp,k2z,k2z,C2);
phase_t = reshape(data(1,1,:,2),[length(omega) 1]);
amp_t = reshape(data(1,1,:,1),[length(omega) 1])*power; % in unit of K
if fit_type == "both"
    y1 = amp_t - amp/CTR;
    y2 = phase_t - phase;
    y = [y1,y2];
elseif fit_type == "phase"
    y= phase_t - phase;
else
    y = amp_t - amp/CTR;
end
end