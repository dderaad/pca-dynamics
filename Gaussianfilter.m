function filter = Gaussianfilter(t,a,tau)

filter = 1/(sqrt(2*pi)*a)*exp(1/(2*a^2)*-(t-tau).^2);

end

