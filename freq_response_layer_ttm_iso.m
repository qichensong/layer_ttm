function Itgd = freq_response_layer_ttm_iso(omega,w0,w1,k,rho_cp,kzz,kxx,alpha,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
	global N_imag_layer
    M = length(k);
	Itgd = zeros(length(omega),M);			% term depending on heat transfer problem in integrand (-D/C)

    L1 = L(1);
    L2 = L(2);
    
    for mm = 1:M
        q = k(mm);
        for iw = 1:length(omega)
            lambdap = (kxx*(q^2)+1i*rho_cp*omega(iw)+g)/kzz;
            lambdae = (kex*(q^2)+1i*Ce*omega(iw)+g)/kez;
            Mmat = [lambdae,-g/kez;-g/kzz,lambdap];
            [u,D] = eig(Mmat);
            lambda1 = D(1,1);
            lambda2 = D(2,2);
            u11 = u(1,1);
            u12 = u(1,2);
            u21 = u(2,1);
            u22 = u(2,2);
            U = u11*u22-u12*u21;
            %v = inv(u);
            v = 1/U ... 
                *[u22 -u12;
                -u21  u11];

            
            v11 = v(1,1);
            v12 = v(1,2);
            v21 = v(2,1);
            v22 = v(2,2);
            
            D1 = v11/(lambda1-alpha^2)/kez;
            D2 = v21/(lambda2-alpha^2)/kez;

            % thickness of the first layer
            zz = L1;
            

            sl1 = sqrt(lambda1);
            sl2 = sqrt(lambda2);
             
% 
%             Q11 = 1/U*(u11*u22*cosh(sl1*z)-u12*u21*cosh(sl2*z));
%             Q12 = u11*u12/U*(-cosh(sl1*z)+cosh(sl2*z));
%             Q13 = u11*(alpha/sl1*sinh(sl1*z)-cosh(sl1*z)+exp(-alpha*z));
%             Q14 = u12*(alpha/sl2*sinh(sl2*z)-cosh(sl2*z)+exp(-alpha*z));
%             Q21 = u21*u22/U*(cosh(sl1*z)-cosh(sl2*z));
%             Q22 = 1/U*(-u12*u21*cosh(sl1*z)+u11*u22*cosh(sl2*z));
%             Q23 = u21*(alpha/sl1*sinh(sl1*z)-cosh(sl1*z)+exp(-alpha*z));
%             Q24 = u22*(alpha/sl2*sinh(sl2*z)-cosh(sl2*z)+exp(-alpha*z));
%             Q31 = kez/U*(-u11*u22*sl1*sinh(sl1*z)+u12*u21*sl2*sinh(sl2*z));
%             Q32 = kez*u11*u12/U*(sl1*sinh(sl1*z)-sl2*sinh(sl2*z));
%             Q33 = kez*u11*(sl1*sinh(sl1*z)-alpha*cosh(sl1*z)+alpha*exp(-alpha*z));
%             Q34 = kez*u12*(sl2*sinh(sl2*z)-alpha*cosh(sl2*z)+alpha*exp(-alpha*z));
%             Q41 = u21*u22*kzz/U*(-sl1*sinh(sl1*z)+sl2*sinh(sl2*z));
%             Q42 = kzz/U*(u12*u21*sl1*sinh(sl1*z)-u11*u22*sl2*sinh(sl2*z));
%             Q43 = u21*kzz*(sl1*sinh(sl1*z)-alpha*cosh(sl1*z)+alpha*exp(-alpha*z));
%             Q44 = u22*kzz*(sl2*sinh(sl2*z)-alpha*cosh(sl2*z)+alpha*exp(-alpha*z));
%             if abs(sl2*z)<50
%             if abs(sl1)>abs(sl2)
%                 scalef = 1/abs(sl2)/kez;
%                             
%                 Q11 = 1/U*(u11*u22-u12*u21*cosh(sl2*z)/cosh(sl1*z));
%                 Q12 = u11*u12/U*(-1+cosh(sl2*z)/cosh(sl1*z));
%                 Q13 = u11*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z)/cosh(sl1*z));
%                 Q14 = u12*(alpha/sl2*sinh(sl2*z)/cosh(sl1*z)-cosh(sl2*z)/cosh(sl1*z)+exp(-alpha*z)/cosh(sl1*z));
%                 Q21 = u21*u22/U*(1-cosh(sl2*z)/cosh(sl1*z));
%                 Q22 = 1/U*(-u12*u21*1+u11*u22*cosh(sl2*z)/cosh(sl1*z));
%                 Q23 = u21*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z)/cosh(sl1*z));
%                 Q24 = u22*(alpha/sl2*sinh(sl2*z)/cosh(sl1*z)-cosh(sl2*z)/cosh(sl1*z)+exp(-alpha*z)/cosh(sl1*z));
%                 Q31 = kez/U*(-u11*u22*sl1*tanh(sl1*z)+u12*u21*sl2*sinh(sl2*z)/cosh(sl1*z));
%                 Q32 = kez*u11*u12/U*(sl1*tanh(sl1*z)-sl2*sinh(sl2*z)/cosh(sl1*z));
%                 Q33 = kez*u11*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z)/cosh(sl1*z));
%                 Q34 = kez*u12*(sl2*sinh(sl2*z)/cosh(sl1*z)-alpha*cosh(sl2*z)/cosh(sl1*z)+alpha*exp(-alpha*z)/cosh(sl1*z));
%                 Q41 = u21*u22*kzz/U*(-sl1*tanh(sl1*z)+sl2*sinh(sl2*z)/cosh(sl1*z));
%                 Q42 = kzz/U*(u12*u21*sl1*tanh(sl1*z)-u11*u22*sl2*sinh(sl2*z)/cosh(sl1*z));
%                 Q43 = u21*kzz*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z))/cosh(sl1*z);
%                 Q44 = u22*kzz*(sl2*sinh(sl2*z)/cosh(sl1*z)-alpha*cosh(sl2*z)/cosh(sl1*z)+alpha*exp(-alpha*z)/cosh(sl1*z));
%             else
                scalef = 1/abs(sl2)/kzz;
                z =  zz/N_imag_layer;
                
                Q11 = 1/U*(u11*u22*cosh(sl1*z)/cosh(sl2*z)-u12*u21);
                Q12 = u11*u12/U*(-cosh(sl1*z)/cosh(sl2*z)+1);
                Q13 = u11*(alpha/sl1*sinh(sl1*z)/cosh(sl2*z)-cosh(sl1*z)/cosh(sl2*z)+exp(-alpha*z)/cosh(sl2*z));
                Q14 = u12*(alpha/sl2*tanh(sl2*z)-1+exp(-alpha*z)/cosh(sl2*z));
                Q21 = u21*u22/U*(cosh(sl1*z)/cosh(sl2*z)-1);
                Q22 = 1/U*(-u12*u21*cosh(sl1*z)/cosh(sl2*z)+u11*u22);
                Q23 = u21*(alpha/sl1*sinh(sl1*z)/cosh(sl2*z)-cosh(sl1*z)/cosh(sl2*z)+exp(-alpha*z)/cosh(sl2*z));
                Q24 = u22*(alpha/sl2*tanh(sl2*z)-1+exp(-alpha*z)/cosh(sl2*z));
                Q31 = kez/U*(-u11*u22*sl1*sinh(sl1*z)/cosh(sl2*z)+u12*u21*sl2*tanh(sl2*z));
                Q32 = kez*u11*u12/U*(sl1*sinh(sl1*z)/cosh(sl2*z)-sl2*tanh(sl2*z));
                Q33 = kez*u11*(sl1*sinh(sl1*z)/cosh(sl2*z)-alpha*cosh(sl1*z)/cosh(sl2*z)+alpha*exp(-alpha*z)/cosh(sl2*z));
                Q34 = kez*u12*(sl2*tanh(sl2*z)-alpha+alpha*exp(-alpha*z)/cosh(sl2*z));
                Q41 = u21*u22*kzz/U*(-sl1*sinh(sl1*z)/cosh(sl2*z)+sl2*tanh(sl2*z));
                Q42 = kzz/U*(u12*u21*sl1*sinh(sl1*z)/cosh(sl2*z)-u11*u22*sl2*tanh(sl2*z));
                Q43 = u21*kzz*(sl1*sinh(sl1*z)/cosh(sl2*z)-alpha*cosh(sl1*z)/cosh(sl2*z)+alpha*exp(-alpha*z)/cosh(sl2*z));
                Q44 = u22*kzz*(sl2*tanh(sl2*z)-alpha+alpha*exp(-alpha*z)/cosh(sl2*z));

                Q = [Q11 Q12 Q13 Q14
                     Q21 Q22 Q23 Q24
                     Q31 Q32 Q33 Q34
                     Q41 Q42 Q43 Q44];
                Q(3,:) = Q(3,:)*scalef;
                Q(4,:) = Q(4,:)*scalef;
                
                zhalf = 0;
                Qhalf0 = [u11*exp(-sl1*zhalf) u11*exp(sl1*zhalf) u12*exp(-sl2*zhalf) u12*exp(sl2*zhalf)
                         u21*exp(-sl1*zhalf) u21*exp(sl1*zhalf) u22*exp(-sl2*zhalf) u22*exp(sl2*zhalf)
                         u11*kez*sl1*exp(-sl1*zhalf) -u11*kez*sl1*exp(sl1*zhalf) u12*kez*sl2*exp(-sl2*zhalf) -u12*kez*sl2*exp(sl2*zhalf)
                         u21*kzz*sl1*exp(-sl1*zhalf) -u21*kzz*sl1*exp(sl1*zhalf) u22*kzz*sl2*exp(-sl2*zhalf) -u22*kzz*sl2*exp(sl2*zhalf)];
                Qhalf0(3,:) = Qhalf0(3,:)*scalef;
                Qhalf0(4,:) = Qhalf0(4,:)*scalef;
                Qhalf_i = inv(Qhalf0);
                
            for iz = 2:N_imag_layer
            
                zhalf = z/N_imag_layer;
                Qhalf = [u11*exp(-sl1*zhalf) u11*exp(sl1*zhalf) u12*exp(-sl2*zhalf) u12*exp(sl2*zhalf)
                         u21*exp(-sl1*zhalf) u21*exp(sl1*zhalf) u22*exp(-sl2*zhalf) u22*exp(sl2*zhalf)
                         u11*kez*sl1*exp(-sl1*zhalf) -u11*kez*sl1*exp(sl1*zhalf) u12*kez*sl2*exp(-sl2*zhalf) -u12*kez*sl2*exp(sl2*zhalf)
                         u21*kzz*sl1*exp(-sl1*zhalf) -u21*kzz*sl1*exp(sl1*zhalf) u22*kzz*sl2*exp(-sl2*zhalf) -u22*kzz*sl2*exp(sl2*zhalf)];

                Qhalf(3,:) = Qhalf(3,:)*scalef;
                Qhalf(4,:) = Qhalf(4,:)*scalef;     
              

                Q = Qhalf*Qhalf_i*Q;
            end
            
            G = [Gep/(Gep+Gpp) Gpp/(Gep+Gpp) -1/(Gep+Gpp)/scalef -1/(Gpp+Gpp)/scalef
                  0 0  1     1];

            GQ = G*Q;
            lambda = (k2x*q^2+1j*C2*omega(iw))/k2z;
            
            scalef1 = sqrt(scalef)/sqrt(k2z*abs(sqrt(lambda)));
                       
            F = [ 1 -tanh(sqrt(lambda)*L2)/k2z/sqrt(lambda)/scalef
                 -k2z*sqrt(lambda)*tanh(sqrt(lambda)*L2)*scalef1 1/scalef*scalef1];

            T = F*GQ;
%             if any(isnan(Q(:))) 
%                 disp(Q)
%             end

%             scale3 = ;
            
            X = scalef^(-1).*Q(3,:)-Gep.*Q(1,:)+Gep.*GQ(1,:);
            %X = scalef^(-1).*Q(4,:)-Gpp.*Q(1,:)+Gpp.*GQ(1,:);

                 
            Itgd(iw,mm) = -((X(1)*T(2,3)-X(3)*T(2,1))*D1+(X(1)*T(2,4)-X(4)*T(2,1))*D2)/(X(1)*T(2,2)-X(2)*T(2,1))...
                     *exp(-q.^2*w0^2/4)/(2*pi)^2;

             
        end
    end
%     for i = 1:length(omega)
%         hold on
%      plot(k,real((Itgd(i,:))),'-')
%      plot(k,imag((Itgd(i,:))),'--')
%     end

end