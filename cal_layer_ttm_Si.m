close all
T = 295;
C = getC('Au',T);
kz = 2;%phonon thermal conductivity of metal part

global N_imag_layer Isotropic
N_imag_layer = 1; % increase this number if you find numerical instability

wmin = 5e2;
wmax = 5e7;
nw = 100;
omega = exp(linspace(log(wmin),log(wmax),nw))*2*pi;

w0 = 3.14e-6;%radius
w1 = 3.14e-6;

CTR = 3.54e-4;
alp = 1/(17e-9);
gamma = 1/(17e-9);
TLAbs = 0.610;

% For Si substrate

Ce = 1e4;%electron volumetric heat capacity [J/m3K]
C = C - Ce;%phonon volumetric heat capacity [J/m3K]
kez = 195 - kz;%electron thermal conductivity of metal part
g = 2.2e16;%volumetric heating rate of electrons on the same side[W/m3K]
L = [116e-9 1e-3];
Isotropic = [1 1]; % istropic 1, anisotropic 0;
fname = '/Users/60.moon/Desktop/Research/Thermal conductivity/FDTR/100522_Au films with Qichen/Si_2nmTi_25_100nmAu_50_Measurement_1_3.14um_17.29mW_21.88C_Sweep_2.txt';

Gep = 0;%thermal conductance of electrons on the metal side coupled with phonons on the semiconductor side
Gpp = 1e8;%phonon-phonon conductance
k2x = 130;%substrate / in-plane thermal conductivity
k2z = k2x;%second layer - substate / cross-place thermal conductivity
C2 = getC(['Si'],T);


data1 = get_data(fname);% load experimental curve
power = get_power(fname)*TLAbs;%[W]

fit_type = "both";% "both" = fitting phase and amplitude, "phase" = fitting only phase
param = [w0,C,kz,kz,alp,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR];
param_fit = [10,13,17];
var_fit = [10];
var_p = 0.20;
x0 = [116e-9,1e8,3.54e-4];
xu = 10*x0;
xl = 0.1*x0;
options = optimset('TolFun',1e-3,'Display','off','TolX',1e-3);
fh=  @(x) dif_te(x,param_fit,data1(:,1)*2*pi,data1(:,2),data1(:,3),param,power,fit_type);
xx = lsqnonlin(fh,x0,xl,xu,options);
for ix = 1:length(xx)
    disp(xx(ix))
end

param1 = param;
for i = 1:length(param_fit)
    param1(param_fit(i)) = xx(i); % update using the best fits
end
[w0,C,kx,kz,alp,Ce,kex,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR]=unpack_param(param1);
param2 = param1;
param2(var_fit) = param2(var_fit)*(1-var_p); 
param3 = param1;
param3(var_fit) = param3(var_fit)*(1+var_p);


data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2z,k2z,C2);
amp = reshape(data(1,1,:,1),[length(omega) 1]);
phase = reshape(data(1,1,:,2),[length(omega) 1]);
data0 = return_data_real_space(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2z,k2z,C2);



[w0,C,kx,kz,alp,Ce,kex,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR]=unpack_param(param2);
data_lower = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2z,k2z,C2);
amp_lower = reshape(data_lower(1,1,:,1),[length(omega) 1]);
phase_lower = reshape(data_lower(1,1,:,2),[length(omega) 1]);
[w0,C,kx,kz,alp,Ce,kex,kez,g,L,Gep,Gpp,k2x,k2z,C2,CTR]=unpack_param(param3);
data_higher = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2z,k2z,C2);
amp_higher = reshape(data_higher(1,1,:,1),[length(omega) 1]);
phase_higher = reshape(data_higher(1,1,:,2),[length(omega) 1]);


figure(1)
hold on
plot(data1(:,1),data1(:,2)/CTR,'o','MarkerSize',15)
plot(omega/2/pi,data0(:,1)*power,'LineStyle','--','DisplayName','w/o e-p','LineWidth',1.5)
plot(omega/2/pi,amp*power,'LineStyle','-','DisplayName','Averaged Data','LineWidth',1.5)
plot(omega/2/pi,amp_lower*power,'LineStyle','-.','DisplayName','Lower Bound','LineWidth',1.5)
plot(omega/2/pi,amp_higher*power,'LineStyle','-.','DisplayName','Upper Bound','LineWidth',1.5)
set(gca, 'XScale', 'log','FontSize',14)
xlabel('Frequency (Hz)','FontSize',17)
ylabel('Amplitude (K)','FontSize',17)
figure(2)
hold on
plot(data1(:,1),data1(:,3),'o','MarkerSize',15)
plot(omega/2/pi,data0(:,2),'LineStyle','-','DisplayName','w/o e-p','LineWidth',1.5)
plot(omega/2/pi,phase,'LineStyle','--','DisplayName','Averaged Data','LineWidth',1.5)
plot(omega/2/pi,phase_lower,'LineStyle','-.','DisplayName','Lower Bound','LineWidth',1.5)
plot(omega/2/pi,phase_higher,'LineStyle','-.','DisplayName','Upper Bound','LineWidth',1.5)
set(gca, 'XScale', 'log','FontSize',14)
xlabel('Frequency (Hz)','FontSize',17)
ylabel('Phase (deg)','FontSize',17)

% g = 7e16;
% data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
% amp = reshape(data(1,1,:,1),[length(omega) 1]);
% figure(1)
% plot(omega/2/pi,amp/amp(1),'LineStyle','--','DisplayName','no e-p','LineWidth',3)
% figure(2)
% plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','--','DisplayName','no e-p','LineWidth',3)
% legend show


% tdelay = linspace(500e-12,8e-9,100);
% datat = return_tdtr_real_space_layer_ttm(tdelay,6e6*2*pi,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
% figure(3)
% amp = reshape(datat(1,1,:,1),[length(tdelay) 1]);
% hold on
% plot(tdelay,amp,'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Delay time (s)')
% ylabel('Amplitude')
% figure(4)
% hold on
% plot(tdelay,reshape(datat(1,1,:,2),[length(tdelay) 1]),'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Delay time (s)')
% ylabel('Phase (deg)')
% 
% 
% 
% 
% 
% Gep1 = 3e8;
% 
% 
% data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
% amp = reshape(data(1,1,:,1),[length(omega) 1]);
% figure(1)
% hold on
% plot(omega/2/pi,amp,'LineStyle','-','DisplayName','with G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Frequency (Hz)')
% figure(2)
% hold on
% plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','-','DisplayName','with G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Frequency (Hz)')
% 
% tdelay = linspace(500e-12,8e-9,100);
% datat = return_tdtr_real_space_layer_ttm(tdelay,6e6*2*pi,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
% figure(3)
% amp = reshape(datat(1,1,:,1),[length(tdelay) 1]);
% hold on
% plot(tdelay,amp,'LineStyle','-','DisplayName','with w/o G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Delay time (s)')
% figure(4)
% hold on
% plot(tdelay,reshape(datat(1,1,:,2),[length(tdelay) 1]),'LineStyle','-','DisplayName','with w/o G_{ep}','LineWidth',1.5)
% set(gca, 'XScale', 'log')
% xlabel('Delay time (s)')
% for i = 1:4
%     figure(i)
%     box on
%     legend boxoff
%     %saveas(gcf,[num2str(i) '.svg'])
% end
% 
% 
% zz = linspace(0,100e-9,300);
% zz2 = linspace(0,100e-9,200);
% [te,tp,t2]=return_data_real_space_z_layer_ttm(omega,w0,w1,[0],[0],zz,zz2,C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
% figure(5)
% iw = 20;
% % plot(zz,te(:,iw,1).*sin(te(:,iw,2)/180*pi))
% hold on
% % plot(zz,tp(:,iw,1).*sin(tp(:,iw,2)/180*pi))
% % plot(zz(end)+zz2,t2(:,iw,1).*sin(t2(:,iw,2)/180*pi))
% plot(zz,te(:,iw,1).*cos(te(:,iw,2)/180*pi),'r--','LineWidth',2)
% plot(zz,tp(:,iw,1).*cos(tp(:,iw,2)/180*pi),'g--','LineWidth',2)
% plot(zz(end)+zz2,t2(:,iw,1).*cos(t2(:,iw,2)/180*pi),'b--','LineWidth',2)
% figure(6)
% plot(zz,te(:,iw,1),'r--','LineWidth',2)
% hold on
% plot(zz,tp(:,iw,1),'g--','LineWidth',2)
% plot(zz(end)+zz2,t2(:,iw,1),'b--','LineWidth',2)
% [te,tp,t2]=return_data_real_space_z_layer_ttm(omega,w0,w1,[0],[0],zz,zz2,C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
% figure(5)
% % plot(zz,te(:,iw,1).*sin(te(:,iw,2)/180*pi))
% hold on
% % plot(zz,tp(:,iw,1).*sin(tp(:,iw,2)/180*pi))
% % plot(zz(end)+zz2,t2(:,iw,1).*sin(t2(:,iw,2)/180*pi))
% plot(zz,te(:,iw,1).*cos(te(:,iw,2)/180*pi),'r-','LineWidth',2)
% plot(zz,tp(:,iw,1).*cos(tp(:,iw,2)/180*pi),'g-','LineWidth',2)
% plot(zz(end)+zz2,t2(:,iw,1).*cos(t2(:,iw,2)/180*pi),'b-','LineWidth',2)
% figure(6)
% plot(zz,te(:,iw,1),'r-','LineWidth',2)
% hold on
% plot(zz,tp(:,iw,1),'g-','LineWidth',2)
% plot(zz(end)+zz2,t2(:,iw,1),'b-','LineWidth',2)