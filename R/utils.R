
array_operator = function(array1, array2, operation){
  if(is.array(array1) == FALSE){
    array1 = array(array1, dim=dim(array2))
  } else if(is.array(array2) == FALSE){
    array2 = array(array2, dim=dim(array1))
  }
  dim1 = dim(array1)
  dim2 = dim(array2)
  if(length(dim1) < length(dim2)){
    array1 = kronecker(array1, array(1, dim = array(1, dim=length(dim2))))
    dim1 = dim(array1)
  } else if(length(dim1) > length(dim2)){
    array2 = kronecker(array2, array(1, dim = array(1, dim=length(dim1))))
    dim2 = dim(array2)
  }
  
  if(min(dim1-dim2) < 0){
    array1_reshaped = kronecker(array1, array(1, dim=dim2-dim1+1))
    array2_reshaped = array2
  } else if(max(dim1-dim2) > 0){
    array1_reshaped = array1
    array2_reshaped = kronecker(array2, array(1, dim=dim1-dim2+1))
  } else{
    array1_reshaped = array1
    array2_reshaped = array2
  }
  if(operation == '*'){
    array_output = array1_reshaped*array2_reshaped
  } else if(operation == '/'){
    array_output = array1_reshaped/array2_reshaped
  } else if(operation == '+'){
    array_output = array1_reshaped+array2_reshaped
  } else if(operation == '-'){
    array_output = array1_reshaped-array2_reshaped
  } else if(operation == '^'){
    array_output = array1_reshaped^array2_reshaped
  }
  return(array_output)
}

sumDims = function(array, dimension){
  dim1 = length(dim(array))
  keep_dims = setdiff(seq(1, dim1, length.out = dim1), c(dimension))
  array_output = apply(array, MARGIN = keep_dims, FUN = 'sum')
  array_output = kronecker(array_output, array(1, dim=array(1, dim=dim1)))
  dim_reshape = seq(1, dim1, length.out = dim1)
  dim_reshape[(dimension+1):length(dim_reshape)] = dim_reshape[dimension:(length(dim_reshape)-1)]
  dim_reshape[dimension] = dim1
  array_output = aperm(array_output, dim_reshape);
  
  return(array_output)
}

labor_productivity = function(param,data,settings){
  # Settings
  outerdiff = Inf
  w = settings$w_init
  iter = 0
  zeta_init = settings$zeta
  nu_init = 0.005
  nu = nu_init
  
  while(outerdiff>tol & iter<10){
    # 1) Labor supply
    # Indirect utility
    w_tr = aperm(array(w, dim=c(875,3,1)), c(3,1,2))
    reptau = kronecker(data$tau, array(1, dim=c(1, 1, param$S)))
    
    w_tr_reptau = array_operator(w_tr^theta1, reptau^(-theta1), '*')
    lambda_ijs_is = array_operator(w_tr_reptau, sumDims(w_tr_reptau, 2), '/')
    lambda_is_i_p = aperm(kronecker(data$lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));
    
    l_prod = array_operator(data$lambda_i, array_operator(lambda_is_i_p, lambda_ijs_is, '*'), '*')
    L_j_p = sumDims(l_prod, 1)
    L_j = aperm(L_j_p, c(2, 3, 1));
    
    z_L = array_operator(data$L_j_data, L_j, '-')
    w = array_operator(w, 1+nu*(z_L/L_j), '*');
    w = w/w[1,1,1];
    outerdiff = max(abs(z_L))
    
    iter = iter+1;
    
    if(iter>100){
      nu = nu_init*10;
    } else if(iter>1000){
      nu = nu_init*50;
    } else if(iter>10000){
      nu = nu_init*100;
    } else if(iter>1000000){
      break
    }
    print(paste(outerdiff, iter))
  }
  return(list(w=w, w_tr=w_tr, lambda_ijs_is=lambda_ijs_is))
}

amenities = function(param,data,settings){
  B = B_init;
  outerdiff = 10000000;
  iter = 0;
  
  w_tr = param$w_tr
  reptau = param$reptau
  
  while(outerdiff>1e-07){
    # Amenities
    W_is_aux_1 = sumDims(array_operator(w_tr^theta1, reptau^(-theta1), '*'), 2)
    W_is_aux_1 = W_is_aux_1^(1/theta1)
    W_is = aperm(W_is_aux_1, c(1,3,2))
    
    B_W_is = array_operator(B, W_is, '*')^(kappa1)
    
    lambda_is_i_m = array_operator(B_W_is, sumDims(B_W_is, 2), '/');
    lambda_is_i_m_B = array_operator(lambda_is_i_m, B, '/'); 
    
    B_upd = array_operator(lambda_is_i, lambda_is_i_m_B, '/');
    B_upd = array_operator(B_upd, B_upd[,1,], '/');
    
    z_B = array_operator(lambda_is_i, lambda_is_i_m, '-');
    diff_B = max(abs(z_B));
    B = 0.2*B_upd+0.8*kronecker(B, array(1, dim=c(1,1,1)));
    B = array_operator(B, B[,1,1], '/');
    outerdiff = diff_B
    
    iter = iter+1;
    
    if(iter>200){
      break
    }
    print(outerdiff)
  }
  return(list(W_is=W_is, B=B))
}

av_income = function(lambda_is_i, B, W_is, kappa1, y_bar, lambda_ijs_is, w, lambda_i){
  lambda_is_i_p = aperm(kronecker(lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));
  y_bar = (sumDims(array_operator(B, W_is, '*')^(kappa1), 2))^(1/kappa1);
  L_j_eff_p = sumDims(array_operator(y_bar, array_operator(lambda_i, array_operator(lambda_is_i_p, array_operator(lambda_ijs_is, L_bar, '*'), '*'), '*'), '*'),1);
  L_j_eff_w = aperm(L_j_eff_p, c(2, 3, 1));
  L_j_eff = array_operator(L_j_eff_w, w, '/');
  return(list(lambda_is_i_p=lambda_is_i_p, L_j_eff_p=L_j_eff_p, L_j_eff_w=L_j_eff_w, L_j_eff=L_j_eff, y_bar=y_bar))
}

floorspace_supply = function(delta1, H_bar_rest, H_barr, w, L_j_eff, y_bar, lambda_i, L_bar, param){
  gamma1=array_operator(delta1, (1-delta1), '/');
  H_r = pmin(H_bar_rest,H_bar);
  H_f = pmin(H_bar_rest,H_bar);
  beta_tilde_sector = aperm(array_operator((1-param$beta1), param$beta1, '/'), c(1, 3, 2));
  q = (array_operator(sumDims(array_operator(beta_tilde_sector, array_operator(w, L_j_eff, '*'), '*'),2), H_f, '/'));
  r = ((1-alpha1)/alpha1)*(array_operator(array_operator(y_bar, lambda_i*L_bar, '*')+array_operator(q, H_f, '*'), H_r, '/'));        
  return(list(beta_tilde_sector=beta_tilde_sector, q=q, r=r, H_r=H_r, H_f=H_f))
}

number_firms = function(beta1, beta_tilde_sector, w, L_j_eff, sigma1, q){
  beta_tilde = aperm(array_operator(array_operator(beta1, (-beta1), '^'), (array_operator((1-beta1), (-(1-beta1)), '^')), '*'), c(1, 3, 2));
  H_f_sector = array_operator((array_operator(kronecker(beta_tilde_sector, array(1, dim=c(875,1,1))), sumDims(array_operator(w, L_j_eff, '*'), 2), '*')), q, '/');
  sigma_cons = aperm(sigma1, c(1, 3, 2)); 
  beta_factor = aperm(beta1, c(1, 3, 2));
  M = array_operator(array_operator(array_operator(beta_tilde, array_operator(L_j_eff, (beta_factor), '^'), '*'), array_operator(H_f_sector,(1-beta_factor), '^'), '*'), array_operator(sigma_cons, F, '*'), '/');        
  return(list(beta_factor=beta_factor, M=M, sigma_cons=sigma_cons))
}

price_location = function(sigma1, beta_factor, w, q, A){
  markup = aperm(array_operator(sigma1, (sigma1-1), '/'), c(1, 3, 2));
  inv_sigma = aperm((1/(1-sigma1)), c(1, 3, 2));
  p_j = array_operator(array_operator(array_operator(markup, array_operator(w, beta_factor, '^'), '*'), array_operator(q, kronecker(1-beta_factor, array(1, dim=c(875,1,1))), '^'), '*'), A, '/');        
  return(list(p_j=p_j, inv_sigma=inv_sigma))
}

agg_price_index = function(sigma1, M, p_j, sigma1_tr, inv_sigma){
  sigma1_tr = aperm(sigma1, c(1, 3, 2));
  P_s = array_operator((sumDims(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'),1)), (inv_sigma), '^');
  pi_js = array_operator(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'), (array_operator(P_s, (1-sigma1_tr), '^')), '/');        
  P = array_operator(sumDims(array_operator(P_s, (1-xi1), '^'),2), (1/(1-xi1)), '^');
  pi_s = array_operator(array_operator(P_s, (1-xi1), '^'), (array_operator(P, (1-xi1), '^')), '/');
  return(list(pi_js=pi_js, pi_s=pi_s, P=P, P_s=P_s))
}

total_sales = function(y_bar, lambda_i, L_bar, r, H_r, q, H_f, alpha1, pi_js, pi_s){
  X = sumDims(array_operator(y_bar, array_operator(lambda_i, kronecker(L_bar, array(1, dim=c(875))), '*'), '*') + array_operator(r, H_r, '*') +array_operator(q, H_f, '*'), 1);
  Y_js = array_operator(alpha1, array_operator(pi_js, array_operator(pi_s, X, '*'), '*'), '*');
  return(list(X=X, Y_js=Y_js))
}

new_labor_demand = function(sigma_cons, w, L_j_eff, beta_factor, Y_js, A, inv_sigma){
  inv_sigma = array_operator(sigma_cons-1, -1, '^');
  LS = array_operator(w, L_j_eff, '*');
  LD = array_operator(beta_factor, Y_js, '*');
  LD_A = array_operator(LD, array_operator(A, (sigma_cons-1), '^'), '/');
  A_prime = array_operator(array_operator(LS, LD_A, '/'), inv_sigma, '^');
  return(list(inv_sigma=inv_sigma, LS=LS, LD=LD, LD_A=LD_A, A_prime=A_prime))
}

eq_quantities = function(param,data,settings){
  outerdiff = 100000000;
  A = A_init;
  iter = 0;
  zeta_init = param$zeta_init;
  zeta = 0.1;
  
  w = param$w
  W_is = param$W_is
  B = param$B
  lambda_ijs_is = param$lambda_ijs_is
  
  # 2) Average income in each location
  av_inc = av_income(lambda_is_i, B, W_is, kappa1, y_bar, lambda_ijs_is, w, lambda_i)
  L_j_eff = av_inc$L_j_eff
  y_bar = av_inc$y_bar
  
  # 3) Floorspace supply
  floorspace_sup = floorspace_supply(delta1, H_bar_rest, H_barr, w, L_j_eff, y_bar, lambda_i, L_bar, param)
  beta_tilde_sector = floorspace_sup$beta_tilde_sector
  q = floorspace_sup$q
  r = floorspace_sup$r
  H_r = floorspace_sup$H_r
  H_f = floorspace_sup$H_f
  
  # 5) Number of firms by location
  number_fir = number_firms(beta1, beta_tilde_sector, w, L_j_eff, sigma1, q)
  beta_factor = number_fir$beta_factor
  M = number_fir$M
  sigma_cons = number_fir$sigma_cons
  
  while(outerdiff>1e-06 & iter < 10){
    # 6) Price in each location
    price_loc = price_location(sigma1, beta_factor, w, q, A)
    p_j = price_loc$p_j
    inv_sigma = price_loc$inv_sigma
    
    # 7) Aggregate price index and expenditure share in the variety of each location
    agg_price_ind = agg_price_index(sigma1, M, p_j, sigma1_tr, inv_sigma)
    pi_js = agg_price_ind$pi_js
    pi_s = agg_price_ind$pi_s
    P = agg_price_ind$P
    
    # 8) Total sales in each location
    total_sal = total_sales(y_bar, lambda_i, L_bar, r, H_r, q, H_f, alpha1, pi_js, pi_s)
    X = total_sal$X
    Y_js = total_sal$Y_js
    
    # 9) New labor demand
    new_labor_dem = new_labor_demand(sigma_cons, w, L_j_eff, beta_factor, Y_js, A, inv_sigma)
    inv_sigma = new_labor_dem$inv_sigma
    LS = new_labor_dem$LS
    LD = new_labor_dem$LD
    LD_A = new_labor_dem$LD_A
    A_prime = new_labor_dem$A_prime
    
    # 10) Update
    z_A = array_operator(A, A_prime, '-');
    A = array_operator(zeta*A_prime, (1-zeta)*A, '+');
    
    A_diff = max(abs(z_A));
    
    outerdiff = max(A_diff)
    iter = iter + 1;
    
    num_U_i = array_operator(sumDims(array_operator(W_is, kappa1, '^'),2), (1/kappa1), '^');
    den_U_i = array_operator(array_operator(r, (1-alpha1), '^'), array_operator(P, alpha1, '^'), '*');
    U_i = array_operator(num_U_i, den_U_i, '/');
    print(outerdiff)
    
    if(outerdiff<500){
      zeta = zeta_init*10;
    }
    if(outerdiff<300){
      zeta = zeta_init*50;
    }
    if(outerdiff<100){
      zeta = zeta_init*100;
    }
    if(outerdiff<20){
      zeta = zeta_init*500;
    }
    if(outerdiff<10){
      zeta = zeta_init*500;
    }
    if(outerdiff<1){
      zeta = zeta_init*500;
    }
    if(iter>100000){
      break
    }
  }
  return(list(U_i, A))
}

last_eq_quantities = function(param,data,settings){
  outerdiff = 100;
  u = u_init;
  iter = 0;
  zeta_init = param$zeta_init;
  zeta = 0.1;
  
  U_i = param$U_i
  
  while(outerdiff>1e-06){
    lambda_i_model = array_operator(array_operator(u, array_operator(U_i, eta1, '^'), '*'), (sumDims(array_operator(u, array_operator(U_i, eta1, '^'), '*'),1)), '/');
    lambda_i_m_u = array_operator(lambda_i_model, u, '/');
    u_prime = array_operator(lambda_i, lambda_i_m_u, '/');
    
    u = array_operator(zeta*u_prime, (1-zeta)*u, '+');
    u = array_operator(u, u[1,1,1], '/');
    z_u = lambda_i_model-lambda_i;
    diff_u = max(abs(z_u));
    outerdiff = diff_u  
    print(outerdiff)
  }
  return(list(u=u))
}

inversionModel_Eff = function(param,data,settings){

  # Labor productivity
  LP = labor_productivity(param,data,settings)
  param$w = LP$w
  param$w_tr = LP$w_tr
  
  param$reptau = kronecker(data$tau, array(1, dim=c(1, 1, param$S)))
  
  # Amenities
  AM = amenities(param,data,settings)
  param$zeta_init = settings$zeta
  param$W_is = AM$W_is
  param$B = AM$B
  param$lambda_ijs_is = LP$lambda_ijs_is
  
  # Other equilibrium quantities
  print("eq quantities starting")
  EQ = eq_quantities(param,data,settings)
  param$U_i = EQ[[1]]
  A = EQ[[2]]
  
  LEQ = last_eq_quantities(param,data,settings)
  u = LEQ$u
  
  # Save and export
  U = (sumDims(array_operator(u, array_operator(param$U_i, eta1, '^'), '*'),1))^(1/eta1)
  # results = v2struct(u,A,w,B,r,q,pi_js,pi_s,lambda_ijs_is,lambda_is_i,H_r,H_f,U,U_i,lambda_i,L_j,LD,LS,y_bar);
  # return(results)
  
  return(list(A=A, u=u, B=param$B, w=param$w))
}








solveModel1_Eff = function(param,data,settings){
  # Solve baseline model

  # Settings
  outerdiff = Inf;
  w = settings$w_init;
  u = data$u_init;
  iter = 0;
  zeta_init = settings$zeta;

  while(outerdiff>tol){
    # 1) Labor supply
    # Indirect utility
    w_tr = aperm(array(w, dim=c(875,3,1)), c(3,1,2))
    reptau = kronecker(data$tau, array(1, dim=c(1, 1, param$S)))
    
    w_tr_reptau = array_operator(w_tr^theta1, reptau^(-theta1), '*')
    lambda_ijs_is = array_operator(w_tr_reptau, sumDims(w_tr_reptau, 2), '/')

    W_is_aux_1 = sumDims(w_tr_reptau,2);
    W_is_aux = W_is_aux_1^(1/theta1);
    W_is = aperm(W_is_aux, c(1,3,2));
    
    B_W_is = array_operator(data$B, W_is, '*')^(kappa1)
    lambda_is_i = array_operator(B_W_is, sumDims(B_W_is,2), '/');
    
    lambda_is_i_p = aperm(kronecker(data$lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));
    
    l_prod = array_operator(data$lambda_i, array_operator(lambda_is_i_p, lambda_ijs_is, '*'), '*')
    L_j_p = sumDims(l_prod, 1)
    L_j = aperm(L_j_p, c(2, 3, 1));

    # 2) Average income in each location
    av_inc = av_income(lambda_is_i, data$B, W_is, kappa1, y_bar, lambda_ijs_is, w, lambda_i);
    lambda_is_i_p=av_inc$lambda_is_i_p
    L_j_eff_p=av_inc$L_j_eff_p
    L_j_eff_w=av_inc$L_j_eff_w
    L_j_eff=av_inc$L_j_eff
    y_bar=av_inc$y_bar

    # 3) Floorspace supply
    floorspace = floorspace_supply(delta1, H_bar_rest, H_barr, w, L_j_eff, y_bar, lambda_i, L_bar, param)
    H_r = floorspace$H_r
    H_f = floorspace$H_f
    beta_tilde_sector = floorspace$beta_tilde_sector
    q = floorspace$q
    r = floorspace$r

    # 5) Number of firms by location
    num_firms = number_firms(beta1, beta_tilde_sector, w, L_j_eff, sigma1, q)
    beta_factor=num_firms$beta_factor
    M=num_firms$M
    sigma_cons=num_firms$sigma_cons
    
    H_f_sector = array_operator(array_operator(kronecker(beta_tilde_sector, array(1, dim=c(875,1,1))), sumDims(array_operator(w, L_j_eff, '*'),2), '*'), q, '/');

    Q_js = array_operator(array_operator(L_j_eff, beta_factor, '^'), array_operator(H_f_sector, 1-beta_factor, '^'), '*');
    elast_demand = aperm(array_operator(sigma1, sigma1-1, '/'), c(1,3,2));
    inv_elast_demand = aperm(array_operator(sigma1-1, sigma1, '/'), c(1,3,2));
    
    Q_s = array_operator(sumDims(array_operator(Q_js, elast_demand, '^'), 1), inv_elast_demand, '^');
    Q = array_operator(sumDims(array_operator(Q_s, (xi1-1)/xi1, '^'),2), xi1/(xi1-1), '^');
    
    # 6) Price in each location
    price_loc = price_location(sigma1, beta_factor, w, q, data$A)
    p_j = price_loc$p_j
    inv_sigma = price_loc$inv_sigma

    # 7) Aggregate price index and expenditure share in the variety of each location
    agg_price_ind = agg_price_index(sigma1, M, p_j, sigma1_tr, inv_sigma)
    pi_js = agg_price_ind$pi_js
    pi_s = agg_price_ind$pi_s
    P = agg_price_ind$P
    P_s = agg_price_ind$P_s
    
    # 8) Total sales in each location
    total_sal = total_sales(y_bar, lambda_i, L_bar, r, H_r, q, H_f, alpha1, pi_js, pi_s)
    X = total_sal$X
    Y_js = total_sal$Y_js

    # 9) New labor demand and wages
    new_labor_dem = new_labor_demand(sigma_cons, w, L_j_eff, beta_factor, Y_js, A, inv_sigma)
    inv_sigma = new_labor_dem$inv_sigma
    LS = new_labor_dem$LS
    LD = new_labor_dem$LD
    w_prime = array_operator(LD, L_j_eff, '/');
    w_prime = array_operator(w_prime, w_prime[1,1,1], '/');

    # 10) Update
    w = zeta*w_prime + (1-zeta)*w;
    w = array_operator(w, w[1,1,1], '/');

    w_diff = max(abs(w - w_prime));

    num_U_i = (sumDims(array_operator(W_is, kappa1, '^'),2))^(1/kappa1);
    den_U_i = array_operator((r^(1-alpha1)), (P^(alpha1)), '*');
    U_i = array_operator(num_U_i, den_U_i, '/');
    U = array_operator(sumDims(array_operator(u, U_i^eta1, '*'),1), 1/eta1, '^');

    if(settings$endo_Lr==1){
      lambda_i_upd = array_operator(array_operator(u, U_i^eta1, '*'), sumDims(array_operator(u, U_i^eta1, '*'),1), '/');
      lambda_diff = max(abs(lambda_i-lambda_i_upd));
      lambda_i = 0.05*lambda_i_upd+(1-0.05)*lambda_i;
    }
    
    outerdiff = max(max(w_diff), max(lambda_diff))
    iter = iter + 1;

    if(outerdiff<500){
      zeta = zeta_init*10;
    }
    
    if(outerdiff<300){
      zeta = zeta_init*50;
    }
    
    if(outerdiff<100){
      zeta = zeta_init*100;
    }
    
    if(outerdiff<20){
      zeta = zeta_init*500;
    }
    
    if(outerdiff<10){
      zeta = zeta_init*500;
    }
    
    if(outerdiff<1){
      zeta = zeta_init*500;
    }
    print(outerdiff)
  }
  
  # Save and export
  M_sum = sumDims(M,1);
  return(w=w, W_is=W_is, B=B, r=r, q=q, pi_js=pi_js, pi_s=pi_s, lambda_ijs_is=lambda_ijs_is,
         y_bar=y_bar, lambda_is_i=lambda_is_i, H_r=H_r, H_f=H_f, U=U, U_i=U_i,
         lambda_i=lambda_i, L_j=L_j, H_bar=H_bar, H_bar_rest=H_bar_rest, Q_js=Q_js,
         Q_s=Q_s, Q=Q, L_j_eff=L_j_eff, M_sum=M_sum)
}
