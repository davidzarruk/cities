
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
  return(list(w, w_tr, lambda_ijs_is))
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
  return(list(W_is, B))
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
  lambda_is_i_p = aperm(kronecker(lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));
  y_bar = (sumDims(array_operator(B, W_is, '*')^(kappa1), 2))^(1/kappa1);
  L_j_eff_p = sumDims(array_operator(y_bar, array_operator(lambda_i, array_operator(lambda_is_i_p, array_operator(lambda_ijs_is, L_bar, '*'), '*'), '*'), '*'),1);
  L_j_eff_w = aperm(L_j_eff_p, c(2, 3, 1));
  L_j_eff = array_operator(L_j_eff_w, w, '/');
  
  # 3) Floorspace supply
  #delta_tilde = (delta1.^(-delta1)).*((1-delta1).^(-(1-delta1)));
  gamma1=array_operator(delta1, (1-delta1), '/');
  H_r = pmin(H_bar_rest,H_bar);
  H_f = pmin(H_bar_rest,H_bar);
  beta_tilde_sector = aperm(array_operator((1-param$beta1), param$beta1, '/'), c(1, 3, 2));
  q = (array_operator(sumDims(array_operator(beta_tilde_sector, array_operator(w, L_j_eff, '*'), '*'),2), H_f, '/'));
  r = ((1-alpha1)/alpha1)*(array_operator(array_operator(y_bar, lambda_i*L_bar, '*')+array_operator(q, H_f, '*'), H_r, '/'));        
  
  # 5) Number of firms by location
  beta_tilde = aperm(array_operator(array_operator(beta1, (-beta1), '^'), (array_operator((1-beta1), (-(1-beta1)), '^')), '*'), c(1, 3, 2));
  H_f_sector = array_operator((array_operator(kronecker(beta_tilde_sector, array(1, dim=c(875,1,1))), sumDims(array_operator(w, L_j_eff, '*'), 2), '*')), q, '/');
  sigma_cons = aperm(sigma1, c(1, 3, 2)); 
  beta_factor = aperm(beta1, c(1, 3, 2));
  M = array_operator(array_operator(array_operator(beta_tilde, array_operator(L_j_eff, (beta_factor), '^'), '*'), array_operator(H_f_sector,(1-beta_factor), '^'), '*'), array_operator(sigma_cons, F, '*'), '/');        
  
  while(outerdiff>1e-06 & iter < 10){
    # 6) Price in each location
    markup = aperm(array_operator(sigma1, (sigma1-1), '/'), c(1, 3, 2));
    inv_sigma = aperm((1/(1-sigma1)), c(1, 3, 2));
    p_j = array_operator(array_operator(array_operator(markup, array_operator(w, beta_factor, '^'), '*'), array_operator(q, kronecker(1-beta_factor, array(1, dim=c(875,1,1))), '^'), '*'), A, '/');        
    
    # 7) Aggregate price index and expenditure share in the variety of each location
    sigma1_tr = aperm(sigma1, c(1, 3, 2));
    P_s = array_operator((sumDims(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'),1)), (inv_sigma), '^');
    pi_js = array_operator(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'), (array_operator(P_s, (1-sigma1_tr), '^')), '/');        
    P = array_operator(sumDims(array_operator(P_s, (1-xi1), '^'),2), (1/(1-xi1)), '^');
    pi_s = array_operator(array_operator(P_s, (1-xi1), '^'), (array_operator(P, (1-xi1), '^')), '/');
    
    # 8) Total sales in each location
    X = sumDims(array_operator(y_bar, array_operator(lambda_i, kronecker(L_bar, array(1, dim=c(875))), '*'), '*') + array_operator(r, H_r, '*') +array_operator(q, H_f, '*'), 1);
    Y_js = array_operator(alpha1, array_operator(pi_js, array_operator(pi_s, X, '*'), '*'), '*');
    
    # 9) New labor demand
    inv_sigma = array_operator(sigma_cons-1, -1, '^');
    LS = array_operator(w, L_j_eff, '*');
    LD = array_operator(beta_factor, Y_js, '*');
    LD_A = array_operator(LD, array_operator(A, (sigma_cons-1), '^'), '/');
    A_prime = array_operator(array_operator(LS, LD_A, '/'), inv_sigma, '^');
    
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
}

inversionModel_Eff = function(param,data,settings){

  # Labor productivity
  LP = labor_productivity(param,data,settings)
  param$w = LP[[1]]
  param$w_tr = LP[[2]]
  
  param$reptau = kronecker(data$tau, array(1, dim=c(1, 1, param$S)))
  
  # Amenities
  AM = amenities(param,data,settings)
  param$zeta_init = settings$zeta
  param$W_is = AM[[1]]
  param$B = AM[[2]]
  param$lambda_ijs_is = LP[[3]]
  
  # Other equilibrium quantities
  EQ = eq_quantities(param,data,settings)
  param$U_i = EQ[[1]]
  A = EQ[[2]]
  
  
  last_eq_quantities(param,data,settings)
  
  # Save and export
  U = (sumDims(array_operator(u, array_operator(param$U_i, eta1, '^'), '*'),1))^(1/eta1)
  # results = v2struct(u,A,w,B,r,q,pi_js,pi_s,lambda_ijs_is,lambda_is_i,H_r,H_f,U,U_i,lambda_i,L_j,LD,LS,y_bar);
  # return(results)
  
  return(list(A, U, param$B, param$w))
}

