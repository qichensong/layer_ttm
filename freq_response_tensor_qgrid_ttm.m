function Itgd = freq_response_tensor_qgrid_ttm(omega,w0,w1,k,rho_cp,kzz,kxx,alpha,gamma,Ce,kez,kex,g)
	M = length(k);
	I = zeros(length(omega),M,M);     %integrand in the function below
	Itgd = zeros(size(I));			% term depending on heat transfer problem in integrand (-D/C)

    
    for m = 1:M
        for mm = 1:M
            qx = k(m);
            qy = k(mm);
            for iw = 1:length(omega)
                lambdap = (kxx*(qx^2+qy^2)+1i*rho_cp*omega(iw)+g)/kzz;
                lambdae = (kex*(qx^2+qy^2)+1i*Ce*omega(iw)+g)/kez;
                Mmat = [lambdae,-g/kez;-g/kzz,lambdap];
                [u,D] = eig(Mmat);
                lambda1 = D(1,1);
                lambda2 = D(2,2);
                %v = inv(u);
                v = 1/(u(1,1)*u(2,2)-u(1,2)*u(2,1)) ... 
                    *[u(2,2) -u(1,2);
                    -u(2,1)  u(1,1)];
                if real(sqrt(lambda1))<0
                    disp(sqrt(lambda1))
                end
                if real(sqrt(lambda2))<0
                    disp(sqrt(lambda2))
                end

                Itgd(iw,m,mm) = alpha*gamma/kez.*((u(1,1).*v(1,1)./(lambda1-alpha^2)+u(1,2).*v(2,1)./(lambda2-alpha^2))/(alpha+gamma) ...
                 -alpha*u(1,1).*v(1,1)/sqrt(lambda1)/(lambda1-alpha^2)/(sqrt(lambda1)+gamma)...
                 -alpha*u(1,2).*v(2,1)/sqrt(lambda2)/(lambda2-alpha^2)/(sqrt(lambda2)+gamma))...
                 *exp(-qx.^2*w0^2/4-qy.^2*w1^2/4)/(2*pi)^2;
            end
        end
    end    
end