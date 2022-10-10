close all
T = 295;
C = getC('Au',T);
kz = 2;%phonon thermal conductivity of metal part



wmin = 5e2;
wmax = 5e7;
nw = 100;
omega = exp(linspace(log(wmin),log(wmax),nw))*2*pi;

w0 = 3.14e-6;%radius
w1 = 3.14e-6;

alp = 1/(16e-9);
gamma = 1/(17e-9);

Ce = 1e4;%electron volumetric heat capacity [J/m3K]
C = C - Ce;%phonon volumetric heat capacity [J/m3K]
kez = 180 - kz;%electron thermal conductivity of metal part
g = 1e16;%volumetric heating rate of electrons on the same side[W/m3K]

L = [120e-9 1e-3];
Gep = 0;%thermal conductance of electrons on the metal side coupled with phonons on the semiconductor side
Gpp = 6e7;%phonon-phonon conductance
k2x = 133;%second layer - substate / in-plane thermal conductivity
k2z = k2x;%substrate / cross-place thermal conductivity
C2 = getC('Si',T);
%g = 5e12;


data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
amp = reshape(data(1,1,:,1),[length(omega) 1]);

data1 = get_data('/Users/60.moon/Desktop/Research/Thermal conductivity/FDTR/100522_Au films with Qichen/Si_2nmTi_25_100nmAu_50_Measurement_1_3.14um_17.29mW_21.88C_Sweep_2.txt');
figure(1)
hold on
plot(data1(:,1),data1(:,2)/data1(1,2),'o')
plot(omega/2/pi,amp/amp(1),'LineStyle','-','DisplayName','with e-p','LineWidth',1.5)
set(gca, 'XScale', 'log','FontSize',14)
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Amplitude','FontSize',18)
figure(2)
hold on
plot(data1(:,1),data1(:,3),'o')
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','-','DisplayName','with e-p','LineWidth',1.5)
set(gca, 'XScale', 'log','FontSize',14)
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Phase (deg)','FontSize',18)

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