function pdf = mvnpdf(x,mu,sigma)
    sizes = size(x); sizes = sizes(1);
    pdf = zeros(sizes, 1);
    for i=1:sizes
      pdf(i) = (1/(sigma*(2*pi)^(1/2))) * exp(-0.5*((x(i)-mu)/sigma)^2);
    end
end