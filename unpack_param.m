function [w0,C,kx,kz,alp,Ce,kex,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR]=unpack_param(param1)
w0 = param1(1);
C = param1(2);
kx = param1(3);
kz = param1(4);
alp = param1(5);
Ce = param1(6);
kex = param1(7);
kez = param1(8);
g = param1(9);
L = [param1(10),param1(11)];
Gep = param1(12);
Gpp = param1(13);
k2x = param1(14);
k2z = param1(15);
C2 = param1(16);
CTR = param1(17);
end