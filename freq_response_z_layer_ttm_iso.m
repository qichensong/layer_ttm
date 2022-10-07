function [tez,tpz,t2z] = freq_response_z_layer_ttm_iso(omega,zz,zz2,w0,w1,k,rho_cp,kzz,kxx,alpha,gamma,Ce,kez,kex,g,L,Gep,Gpp,k2x,k2z,C2)
	M = length(k);
	Itgd = zeros(length(omega),M);			% term depending on heat transfer problem in integrand (-D/C)
    tez = zeros(length(zz),length(omega),M);
    tpz = zeros(length(zz),length(omega),M);
    t2z = zeros(length(zz2),length(omega),M);
    t1 = zeros(length(zz),1);
    t2 = zeros(length(zz),1);
    
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
            z = L1;
            

            sl1 = sqrt(lambda1);
            sl2 = sqrt(lambda2);
            
            scalef = 1/abs(sl1)/kez;
            scalef1 = 1; %scalef;
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

            Q11 = 1/U*(u11*u22-u12*u21*cosh(sl2*z)/cosh(sl1*z));
            Q12 = u11*u12/U*(-1+cosh(sl2*z)/cosh(sl1*z));
            Q13 = u11*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z)/cosh(sl1*z));
            Q14 = u12*(alpha/sl2*sinh(sl2*z)/cosh(sl1*z)-cosh(sl2*z)/cosh(sl1*z)+exp(-alpha*z)/cosh(sl1*z));
            Q21 = u21*u22/U*(1-cosh(sl2*z)/cosh(sl1*z));
            Q22 = 1/U*(-u12*u21*1+u11*u22*cosh(sl2*z)/cosh(sl1*z));
            Q23 = u21*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z)/cosh(sl1*z));
            Q24 = u22*(alpha/sl2*sinh(sl2*z)/cosh(sl1*z)-cosh(sl2*z)/cosh(sl1*z)+exp(-alpha*z)/cosh(sl1*z));
            Q31 = kez/U*(-u11*u22*sl1*tanh(sl1*z)+u12*u21*sl2*sinh(sl2*z)/cosh(sl1*z));
            Q32 = kez*u11*u12/U*(sl1*tanh(sl1*z)-sl2*sinh(sl2*z)/cosh(sl1*z));
            Q33 = kez*u11*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z)/cosh(sl1*z));
            Q34 = kez*u12*(sl2*sinh(sl2*z)/cosh(sl1*z)-alpha*cosh(sl2*z)/cosh(sl1*z)+alpha*exp(-alpha*z)/cosh(sl1*z));
            Q41 = u21*u22*kzz/U*(-sl1*tanh(sl1*z)+sl2*sinh(sl2*z)/cosh(sl1*z));
            Q42 = kzz/U*(u12*u21*sl1*tanh(sl1*z)-u11*u22*sl2*sinh(sl2*z)/cosh(sl1*z));
            Q43 = u21*kzz*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z))/cosh(sl1*z);
            Q44 = u22*kzz*(sl2*sinh(sl2*z)/cosh(sl1*z)-alpha*cosh(sl2*z)/cosh(sl1*z)+alpha*exp(-alpha*z)/cosh(sl1*z));


%%%     cosh(x) -> exp(x)/2 and sinh(x) -> exp(x)/2 when |x| is large.            
%             Q11 = 1/U*(u11*u22-u12*u21*exp(sl2*z-sl1*z));
%             Q12 = u11*u12/U*(-1+exp(sl2*z-sl1*z));
%             Q13 = u11*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z-sl1*z));
%             Q14 = u12*(alpha/sl2*exp(sl2*z-sl1*z)-exp(sl2*z-sl1*z)+exp(-alpha*z-sl1*z));
%             Q21 = u21*u22/U*(1-exp(sl2*z-sl1*z));
%             Q22 = 1/U*(-u12*u21*1+u11*u22*exp(sl2*z-sl1*z));
%             Q23 = u21*(alpha/sl1*tanh(sl1*z)-1+exp(-alpha*z-sl1*z));
%             Q24 = u22*(alpha/sl2*exp(sl2*z-sl1*z)-exp(sl2*z-sl1*z)+exp(-alpha*z-sl1*z));
%             Q31 = kez/U*(-u11*u22*sl1*tanh(sl1*z)+u12*u21*sl2*exp(sl2*z-sl1*z));
%             Q32 = kez*u11*u12/U*(sl1*tanh(sl1*z)-sl2*exp(sl2*z-sl1*z));
%             Q33 = kez*u11*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z-sl1*z));
%             Q34 = kez*u12*(sl2*exp(sl2*z-sl1*z)-alpha*exp(sl2*z-sl1*z)+alpha*exp(-alpha*z-sl1*z));
%             Q41 = u21*u22*kzz/U*(-sl1*tanh(sl1*z)+sl2*exp(sl2*z-sl1*z));
%             Q42 = kzz/U*(u12*u21*sl1*tanh(sl1*z)-u11*u22*sl2*exp(sl2*z-sl1*z));
%             Q43 = u21*kzz*(sl1*tanh(sl1*z)-alpha*1+alpha*exp(-alpha*z-sl1*z));
%             Q44 = u22*kzz*(sl2*exp(sl2*z-sl1*z)-alpha*exp(sl2*z-sl1*z)+alpha*exp(-alpha*z-sl1*z));

            Q = [Q11 Q12 Q13 Q14
                 Q21 Q22 Q23 Q24
                 Q31 Q32 Q33 Q34
                 Q41 Q42 Q43 Q44];
            Q(3,:) = Q(3,:)*scalef;
            Q(4,:) = Q(4,:)*scalef;
     

            G = [Gep/(Gep+Gpp) Gpp/(Gep+Gpp) -1/(Gep+Gpp)/scalef -1/(Gpp+Gpp)/scalef
                  0 0  1     1];
%             G = [1 1 -1/Gep -1/Gpp
%                   0 0  1      1];
            GQ = G*Q;
            lambda = (k2x*q^2+1j*C2*omega(iw))/k2z;
            F = [ 1 -tanh(sqrt(lambda)*L2)/k2z/sqrt(lambda)/scalef
                 -k2z*sqrt(lambda)*tanh(sqrt(lambda)*L2)*scalef1 1/scalef*scalef1];

            T = F*GQ;
            
            X = scalef^(-1).*Q(3,:)-Gep.*Q(1,:)+Gep.*GQ(1,:);

%             Itgd(iw,mm) = alpha*gamma/kez.*((u(2,1).*v(1,1)./(lambda1-alpha^2)+u(2,2).*v(2,1)./(lambda2-alpha^2))/(alpha+gamma) ...
%                  -alpha*u(2,1).*v(1,1)/sqrt(lambda1)/(lambda1-alpha^2)/(sqrt(lambda1)+gamma)...
%                  -alpha*u(2,2).*v(2,1)/sqrt(lambda2)/(lambda2-alpha^2)/(sqrt(lambda2)+gamma))...
%                  *exp(-q.^2*w0^2/4)/(2*pi)^2;
%             Itgd(iw,mm) = -alpha/kez/U*...
%                 (((T(1,1)*T(2,3)-T(1,3)*T(2,1))/(T(1,1)*T(2,2)-T(1,2)*T(2,1)))*(u12*v21/(alpha^2-lambda2)-u22*v11/(alpha^2-lambda1))...
%                 -((T(1,1)*T(2,4)-T(1,4)*T(2,1))/(T(1,1)*T(2,2)-T(1,2)*T(2,1)))*(u11*v21/(alpha^2-lambda2)-u21*v11/(alpha^2-lambda1)))...
%                 *exp(-q.^2*w0^2/4)/(2*pi)^2;


 
%              Itgd(iw,mm) = -1/kez*...
%                      (((T(1,1)*T(2,3)-T(1,3)*T(2,1))*u12*v21-(T(1,1)*T(2,4)-T(1,4)*T(2,1))*u11*v21)/((T(1,1)*T(2,2)-T(1,2)*T(2,1))*(alpha^2-lambda2)*U)...
%                      -((T(1,2)*T(2,3)-T(1,3)*T(2,1))*u22*v11-(T(1,1)*T(2,4)-T(1,4)*T(2,1))*u21*v11)/((T(1,1)*T(2,2)-T(1,2)*T(2,1))*(alpha^2-lambda1)*U))...
%                      *exp(-q.^2*w0^2/4)/(2*pi)^2;


              te0 = ((X(2)*T(2,3)-X(3)*T(2,2))*D1+(X(2)*T(2,4)-X(4)*T(2,2))*D2)/(X(1)*T(2,2)-X(2)*T(2,1));
              tp0 = -((X(1)*T(2,3)-X(3)*T(2,1))*D1+(X(1)*T(2,4)-X(4)*T(2,1))*D2)/(X(1)*T(2,2)-X(2)*T(2,1));
              
              B1 =  u22/2/U*te0-u12/2/U*tp0-0.5*(1+alpha/sl1)*D1;
              B2 = -u21/2/U*te0+u11/2/U*tp0-0.5*(1+alpha/sl2)*D2;
              c1 = B1 + alpha./sl1*D1;
              c2 = B2 + alpha./sl2*D2;

              
              t1(:) = B1*exp(-sl1.*zz)+c1*exp(sl1.*zz)+D1*exp(-alpha.*zz);
              t2(:) = B2*exp(-sl2.*zz)+c2*exp(sl2.*zz)+D2*exp(-alpha.*zz);
              q1 = sl1.*B1-sl1.*c1+alpha.*D1;
              q2 = sl2.*B2-sl2.*c2+alpha.*D2;
              qe0 = kez*(u11*q1+u12*q2);
              qp0 = kzz*(u21*q1+u22*q2);
              tez(:,iw,mm) = (u11*t1+u12*t2)*exp(-q.^2*w0^2/4)/(2*pi)^2;
              tpz(:,iw,mm) = (u21*t1+u22*t2)*exp(-q.^2*w0^2/4)/(2*pi)^2;
              BB =  GQ*[te0;tp0;D1;D2].*cosh(sl1*z);
              B1 = 0.5*(BB(1)+BB(2)/k2z/sqrt(lambda)/scalef);
              c1 = 0.5*(BB(1)-BB(2)/k2z/sqrt(lambda)/scalef);
              t2z(:,iw,mm) = (B1.*exp(-sqrt(lambda).*zz2)+c1.*exp(sqrt(lambda).*zz2))...
                  *exp(-q.^2*w0^2/4)/(2*pi)^2;          
        end
    end
end