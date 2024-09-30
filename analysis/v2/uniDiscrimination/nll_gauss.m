function nLL = nll_gauss(lapse, mu, sigma, D)
       
    pc  = normcdf(D.comparison_loc, D.standard_loc + mu, sigma)*(1-lapse)+lapse/2;
    nLL =  -sum(log(pc.^D.n_r_response .* (1-pc).^D.n_l_response));

 end   
    
    