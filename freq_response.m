function fr = freq_response(omega,w0,w1,sep,rho_cp,Sr,Sz,d,alpha)
    % Frequency response (including radial heat transfer) of a multilayer stack
	% to a periodic heat input on the top layer.
	%
	% Inputs are angular frequency (a vector), pump beam radius (scalar),
	% probe beam radius (scalar), and the following vectors with length equal
	% to the numbers of layers, starting at the top surface: volumetric heat capacity,
	% in-plane thermal conductivity, out-of-plane thermal conductivity, and
	% layer thickness
	%
	% Aaron Schmidt 03-29-2007
    nkps = 50;
	[k,wk]=lgwt(nkps,0,6/sqrt(w0^2));         %Gauss nodes and weights for k
	M = length(k);
	I = zeros(length(omega),M);     %integrand in the function below
	DC = zeros(size(I));			% term depending on heat transfer problem in integrand (-D/C)
    N = length(rho_cp);             %number of layers (including substrate)    
       
	for m = 1:M        
        A = reshape(ones([length(omega) 1]),size(omega));
        B = reshape(zeros([length(omega) 1]),size(omega));
        C = reshape(zeros([length(omega) 1]),size(omega));
        D = reshape(ones([length(omega) 1]),size(omega));
        b = sqrt(k(m)^2*Sr(1)./Sz(1)+1i.*omega*rho_cp(1)/Sz(1));
        b1 = b; % keep it for later 
        % 1 2
        % 3 4
        temp1 = reshape(ones([length(omega) 1]),size(omega));
        temp2 = -tanh(d(1).*b1)/Sz(1)./b1+1/Sz(1)/alpha-exp(-alpha*d(1))/(Sz(1)*alpha)./cosh(d(1).*b1);
        temp3 = -Sz(1).*b1.*tanh(d(1).*b1);
        temp4 = 1-b1/alpha.*tanh(b1*d(1))-exp(-alpha*d(1))./cosh(b1*d(1));
        
        A0 = A;
        B0 = B;
        C0 = C;
        D0 = D;

        A = temp1.* A0 + temp2.* C0;
        B = temp1.* B0 + temp2.* D0;
        C = temp3.* A0 + temp4.* C0;
        D = temp3.* B0 + temp4.* D0;
        for n = 2:N
            if rho_cp(n) ~= 0
                b = sqrt(k(m)^2*Sr(n)./Sz(n)+1i.*omega*rho_cp(n)/Sz(n));
                temp1 = reshape(ones([length(omega) 1]),size(omega));
                temp2 = -tanh(d(n).*b)/Sz(n)./b;
                temp3 = -Sz(n).*b.*tanh(d(n).*b);
                temp4 = reshape(ones([length(omega) 1]),size(omega));
            else
                temp1 = reshape(ones([length(omega) 1]),size(omega));
                temp2 = -reshape(ones([length(omega) 1]),size(omega))/Sz(n);     %contact resistance (R = K/m^2-W)
                temp3 = reshape(zeros([length(omega) 1]),size(omega));
                temp4 = reshape(ones([length(omega) 1]),size(omega));
            end
            A0 = A;
            B0 = B;
            C0 = C;
            D0 = D;    
            
        A = temp1.* A0 + temp2.* C0;
        B = temp1.* B0 + temp2.* D0;
        C = temp3.* A0 + temp4.* C0;
        D = temp3.* B0 + temp4.* D0;
        end
		DC(:,m) = -B./A*alpha./(alpha^2-b1.^2);
	end
    if sep==0
        for m = 1:M
            I(:,m) = k(m)*DC(:,m)*exp( -k(m)^2*(w0^2)/4)*alpha/2/pi; % w0, effective radius u0^2+u1^2 = 2*w0^2
        end
        fr = I*wk;    % do the integration to get back to real space
        fr = reshape(fr,size(omega));
    else
        int_range = 3*w0;	% integration range
		num_nodes = nkps;		% num of nodes. for spots < 2um, may need to be increased
		[x_pts,wx] = lgwt(num_nodes,-int_range,int_range);  % nodes and weights for x
		[y_pts,wy] = lgwt(floor(num_nodes),0,int_range);	% nodes and weights for y
        Tx = zeros(length(omega),length(x_pts));
		Ty = zeros(length(omega),length(y_pts));
		for n = 1:length(x_pts)
			for m = 1:length(y_pts)
				r = sqrt((x_pts(n)-sep)^2+y_pts(m)^2);	% r in the pump frame
				for l = 1:length(k)
					I(:,l) = k(l)*besselj(0,k(l)*r)*DC(:,l)*exp(-k(l)^2*(w0^2)/8); %  effective radius/sqrt(2)
                end
				r = sqrt(x_pts(n)^2+y_pts(m)^2);	% r in the probe frame
				Ty(:,m) = exp(-2*r^2/(w0^2))*I*wk;    % the inner integral 
            end
			Tx(:,n) = 2*Ty*wy; % this is fourier transfrom from -infty to + infty, thus a factor of 2 is needed
		end
		fr = 2*Tx*wx/(w0^2)/pi;
        fr = reshape(fr,size(omega));     % match input shape
    end
end