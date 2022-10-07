function Itgd = freq_response_tensor_qgrid(omega,w0,w1,k,rho_cp,kzz,kxx,alpha,gamma,tau)
	M = length(k);
	I = zeros(length(omega),M,M);     %integrand in the function below
	Itgd = zeros(size(I));			% term depending on heat transfer problem in integrand (-D/C)
    c2z = kzz/rho_cp/tau;
    c2p = kxx/rho_cp/tau;
    
    for m = 1:M
        for mm = 1:M
            qx = k(m);
            qy = k(mm);
            lambda = (c2p*(qx^2+qy^2)-omega.^2+1i*omega/tau)/c2z;
		    Itgd(:,m,mm) = alpha*gamma*(1+1i*omega*tau)/kzz./(lambda-alpha^2)...
             .*(1/(alpha+gamma)-alpha./sqrt(lambda)./(sqrt(lambda)+gamma))*exp(-qx.^2*w0^2/4-qy.^2*w1^2/4)/(2*pi)^2;
        end
    end    
end