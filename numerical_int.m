function fr = numerical_int(irec,omega,k,wk,x,y)
M = length(k);    
for i = 1:size(x,1) %nd
    for ii = 1:size(x,2)
        Tx = zeros(length(omega),M);
        sepx = x(i,ii);
        sepy = y(i,ii);
        for n = 1:M
            Tx(:,n) = reshape(irec(:,n,:),[length(omega),M])*(exp(-1i*k*sepy).*wk);
        end
        frtemp = Tx*(exp(-1i*k*sepx).*wk);
        fr(:,i,ii) = reshape(frtemp,size(omega));     % match input shape
    end
end
end