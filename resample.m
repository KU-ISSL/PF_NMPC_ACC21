function [uk, wk, idx] = resample(uk, wk, resampling_strategy)

Ns = length(wk);  % Ns = number of particles

switch resampling_strategy
   case 'multinomial_resampling'
      with_replacement = true;
      idx = randsample(1:Ns, Ns, with_replacement, wk);
   case 'systematic_resampling'
      edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
      edges(end) = 1;                 % get the upper edge exact
      u1 = rand/Ns;
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
      [~, idx] = histc(u1:1/Ns:1, edges);
   otherwise
      error('Resampling strategy not implemented')
end

uk = uk(:,idx);  
wk = repmat(1/Ns, Ns, 1);          % set all particles have the same weight

end


