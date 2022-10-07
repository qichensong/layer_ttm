function fr = numerical_int_hankel(irec,omega,k,wk)
M = length(k);    
Tx = zeros(length(omega),M);
for n = 1:length(omega)
    Tx(n,:) = irec(n,:).*k(:)';
end
frtemp = Tx*wk*2*pi;
fr(:) = reshape(frtemp,size(omega));     % match input shape
end