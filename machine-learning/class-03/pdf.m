function y = pdf(x, mean, cov)
  y= 1/(2*pi*det(cov)^(1/2))*exp(-0.5.*(x-mean)*inv(cov)*(x-mean)');
end

