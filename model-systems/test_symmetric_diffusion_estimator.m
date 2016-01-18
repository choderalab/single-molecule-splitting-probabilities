for sample = 1:nsamples
  p_i = p_i_samples(sample,:);

  for i = 1:(nbins-1)
    %D_i(i) = dx^2 * Kij(i,i+1) * sqrt(p_i(i) / p_i(i+1));
    D_i(i) = dx^2 * 0.5 * (Kij(i,i+1) * sqrt(p_i(i) / p_i(i+1)) + Kij(i+1,i) * sqrt(p_i(i+1) / p_i(i)));
  end  
  D_i

  % Store samples.
  D_i_samples(sample,:) = D_i;
end

