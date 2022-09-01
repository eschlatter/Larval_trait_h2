model_run <- function(param_list,response,model_statement,n_vars,n_mcmc_iter){
  #structure for plugging in the params to generate a full prior -- this is what MCMCglmm takes
  G <- list()
  for(i in 1:n_vars){
    varname <- paste("G",i,sep="")
    G[[varname]] <- param_list
  }
  if(length(param_list)==2) R=param_list else R=list(V = 1,nu = .002) #set prior for residual
  prior <- list(R=R, #variance prior for residual
                G=G)  #variance priors for random effects
  
  #run the model
  model <- MCMCglmm(formula(paste0(response,'~','1+age')),
                    random=formula(paste0('~',model_statement)),
                    prior=prior, ginverse=list(animal=Ainv),
                    data=phens,nitt=n_mcmc_iter,burnin=round(.01*n_mcmc_iter),thin=10,verbose=FALSE)
  #I think the "missing values in random predictors" warning is just because the parents in phens don't have predictor values
  return(model)
}


#function that takes an mcmc model object (output of MCMCglmm) and looks at the results
model_diagnostics <- function(model){
  #trace and density plots
  plot(model)
  
  #proportional variances as new mcmc objects
  scaled_vars <- list()
  for(i in 1:ncol(model$VCV)){
    varname <- paste('scale',colnames(model$VCV)[i],sep='_')
    scaled_vars[[varname]] <- model$VCV[,i]/rowSums(as.data.frame(model$VCV)) #store them in a list
  }
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  
  #store convergence info by variable (original and scaled)
  diagnostics=data.frame()
  #original variables
  for(i in 1:ncol(model$VCV)){
    v_name <- colnames(model$VCV)[i]
    v_mcmc <- model$VCV[,i]
    v_diagnostics <- data.frame(variable=v_name,
                                effSamp=effectiveSize(v_mcmc),
                                firstAutocorr=autocorr(v_mcmc)[2],lastAutocorr=autocorr(v_mcmc)[nrow(autocorr(v_mcmc))])
    diagnostics=rbind(diagnostics,v_diagnostics)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- scaled_vars[[i]]
    v_diagnostics <- data.frame(variable=v_name,
                                effSamp=effectiveSize(v_mcmc),
                                firstAutocorr=autocorr(v_mcmc)[2],lastAutocorr=autocorr(v_mcmc)[nrow(autocorr(v_mcmc))])
    diagnostics=rbind(diagnostics,v_diagnostics)
  }
  rownames(diagnostics)=c(colnames(model$VCV),colnames(matrix_scaled_vars))
  print(diagnostics)
  return(diagnostics)
}


model_results <- function(model){
  #proportional variances as new mcmc objects
  scaled_vars <- list()
  for(i in 1:ncol(model$VCV)){
    varname <- paste('scale',colnames(model$VCV)[i],sep='_')
    scaled_vars[[varname]] <- model$VCV[,i]/rowSums(as.data.frame(model$VCV)) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  plot(mcmc_areas(matrix_scaled_vars,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:ncol(model$VCV)){
    v_name <- colnames(model$VCV)[i]
    v_mcmc <- model$VCV[,i]
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- scaled_vars[[i]]
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model$VCV),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
  return(estimates)
  
}

#same as above function, but takes a dataframe as input; last column should be chain number
model_results_df <- function(model_df){
  #proportional variances as new mcmc objects
  num_vars <- ncol(model_df)-1
  scaled_vars <- list()
  for(i in 1:num_vars){
    varname <- paste('scale',colnames(model_df[i]),sep='_')
    scaled_vars[[varname]] <- model_df[,i]/rowSums(as.data.frame(model_df[1:num_vars])) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  plot(mcmc_areas(matrix_scaled_vars,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:num_vars){
    v_name <- colnames(model_df)[i]
    v_mcmc <- as.mcmc(model_df[,i])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- as.mcmc(scaled_vars[[i]])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model_df[1:num_vars]),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
  return(estimates)
  
}


#same as above function, but takes a list of 3 model objects as input; should be 3 repeated runs (chains) of the same model
model_results_chains <- function(model_chains){
  chain1 <- as.data.frame(model_chains[[1]]$VCV)
  chain1$chain <- 1
  chain2 <- as.data.frame(model_chains[[2]]$VCV)
  chain2$chain <- 2
  chain3 <- as.data.frame(model_chains[[3]]$VCV)
  chain3$chain <- 3
  
  model_df <- rbind(chain1,chain2,chain3)
  
  #proportional variances as new mcmc objects
  num_vars <- ncol(model_df)-1
  scaled_vars <- list()
  for(i in 1:num_vars){
    varname <- paste('scale',colnames(model_df[i]),sep='_')
    scaled_vars[[varname]] <- model_df[,i]/rowSums(as.data.frame(model_df[1:num_vars])) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  
  #add dummy columns for dam and sire, if necessary
  if(num_vars==3){
    matrix_scaled_vars_plot <- cbind(matrix_scaled_vars,rep(0,nrow(matrix_scaled_vars)),rep(0,nrow(matrix_scaled_vars)))
    matrix_scaled_vars_plot <- matrix_scaled_vars_plot[,c(1,2,4,5,3)]
    colnames(matrix_scaled_vars_plot)[c(3,4)] <- c('scale_dam',"scale_sire")
  }
  else matrix_scaled_vars_plot <- matrix_scaled_vars
  
  plot(mcmc_areas(matrix_scaled_vars_plot,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:num_vars){
    v_name <- colnames(model_df)[i]
    v_mcmc <- as.mcmc(model_df[,i])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- as.mcmc(scaled_vars[[i]])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model_df[1:num_vars]),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
}

#function to run and plot the gelman-rubin diagnostic (using gelman-rubin functions from the coda package)
#input: a list of 3 model objects (should be 3 repeated runs/chains of the same model)
run_gelman_diagnostics <- function(model_chains){
  chain1 <- as.data.frame(model_chains[[1]]$VCV)
  chain1$chain <- 1
  chain2 <- as.data.frame(model_chains[[2]]$VCV)
  chain2$chain <- 2
  chain3 <- as.data.frame(model_chains[[3]]$VCV)
  chain3$chain <- 3
  
  model_chains_df <- rbind(chain1,chain2,chain3)
  mcmc_trace(model_chains_df)
  
  ## Gelman-Rubin diagnostic: looks good; upper limit is basically 1
  
  intercept_gel=mcmc.list(model_chains[[1]]$Sol,model_chains[[2]]$Sol,model_chains[[3]]$Sol)
  clutch_gel=mcmc.list(model_chains[[1]]$VCV[,'clutch'],model_chains[[2]]$VCV[,'clutch'],model_chains[[3]]$VCV[,'clutch'])
  animal_gel=mcmc.list(model_chains[[1]]$VCV[,'animal'],model_chains[[2]]$VCV[,'animal'],model_chains[[3]]$VCV[,'animal'])
  units_gel=mcmc.list(model_chains[[1]]$VCV[,'units'],model_chains[[2]]$VCV[,'units'],model_chains[[3]]$VCV[,'units'])
  if(ncol(model_chains[[1]]$VCV)==5){
    dam_gel=mcmc.list(model_chains[[1]]$VCV[,'dam'],model_chains[[2]]$VCV[,'dam'],model_chains[[3]]$VCV[,'dam'])
    sire_gel=mcmc.list(model_chains[[1]]$VCV[,'sire'],model_chains[[2]]$VCV[,'sire'],model_chains[[3]]$VCV[,'sire'])
  }
  
  #convergence diagnostic values (potential scale reduction factor):
  #"Approximate convergence is diagnosed when the upper limit is close to 1"
  #https://rdrr.io/cran/coda/man/gelman.diag.html
  print(gelman.diag(intercept_gel))
  print(gelman.diag(clutch_gel))
  print(gelman.diag(animal_gel))
  print(gelman.diag(units_gel))
  if(ncol(model_chains[[1]]$VCV)==5){
    print(gelman.diag(dam_gel))
    print(gelman.diag(sire_gel))
  }
  
  #accompanying plots
  gelman.plot(intercept_gel)
  gelman.plot(clutch_gel,main="clutch")
  gelman.plot(animal_gel,main="animal")
  gelman.plot(units_gel,main="units")
  if(ncol(model_chains[[1]]$VCV)==5){
    gelman.plot(dam_gel,main='dam')
    gelman.plot(sire_gel,main='sire')
  }
}


#function to plot prior and posterior
plot_priorandpost <- function(model,nu.,V.,xlim=0.4){
  
  #calculate modes of both distributions
  xvals = seq(0,xlim,length.out=nrow(model$VCV))
  prior_mode = xvals[which.max(dinvgamma(xvals,shape = nu./2,scale = nu.*V./2))]
  
  #plot (VA)
  animal_samples = as.data.frame(model$VCV)$animal
  post_density = density(animal_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  
  colors=c('Prior'='blue','Posterior'='red')
  xlim=max(post_density$x)
  g <- ggplot(as.data.frame(model$VCV), aes(x=animal)) + 
    geom_density(aes(color='Posterior'))+
    geom_line(aes(x = seq(0,xlim,length.out=nrow(model$VCV)), 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V./2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('animal')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g)
  
  #plot (VC)
  clutch_samples = as.data.frame(model$VCV)$clutch
  post_density = density(clutch_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  xlim=max(post_density$x)
  colors=c('Prior'='blue','Posterior'='red')
  
  g <- ggplot(as.data.frame(model$VCV), aes(x=clutch)) + 
    geom_density(aes(color='Posterior'))+
    geom_line(aes(x = seq(0,xlim,length.out=nrow(model$VCV)), 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V./2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('clutch')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g)
  
  #plot (VE)
  resid_samples = as.data.frame(model$VCV)$units
  post_density = density(resid_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  xlim=max(post_density$x)
  colors=c('Prior'='blue','Posterior'='red')
  
  g <- ggplot(as.data.frame(model$VCV), aes(x=units)) + 
    geom_density(aes(color='Posterior'))+
    geom_line(aes(x = seq(0,xlim,length.out=nrow(model$VCV)), 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V./2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('units')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g) 
}

#function to plot prior and posterior
plot_priorandpost_new <- function(model,nu.,V_A,V_C,V_E,xlim=0.4){
  
  xvals = seq(0,xlim,length.out=nrow(model$VCV))  
  
  #VA
  prior_mode = xvals[which.max(dinvgamma(xvals,shape = nu./2,scale = nu.*V_A/2))]
  animal_samples = as.data.frame(model$VCV)$animal
  post_density = density(animal_samples)
  post_density$y = post_density$y/length(animal_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  
  colors=c('Prior'='blue','Posterior'='red')
  xlim=max(post_density$x)
  g <- ggplot(as.data.frame(model$VCV), aes(x=animal)) + 
    geom_density(aes(y=..scaled..,color='Posterior'))+
    geom_line(aes(x = seq(0,xlim,length.out=nrow(model$VCV)), 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V_A/2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('animal')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g)
  
  #VC
  prior_mode = xvals[which.max(dinvgamma(xvals,shape = nu./2,scale = nu.*V_C/2))]
  clutch_samples = as.data.frame(model$VCV)$clutch
  post_density = density(clutch_samples)
  post_density$y = post_density$y/length(clutch_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  xlim=max(post_density$x)
  colors=c('Prior'='blue','Posterior'='red')
  
  g <- ggplot(as.data.frame(model$VCV), aes(x=clutch)) + 
    geom_density(aes(y=..scaled..,color='Posterior'))+
    geom_line(aes(x = seq(0,xlim,length.out=nrow(model$VCV)), 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V_C/2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('clutch')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g)
  
  #plot (VE)
  prior_mode = xvals[which.max(dinvgamma(xvals,shape = nu./2,scale = nu.*V_E/2))]
  resid_samples = as.data.frame(model$VCV)$units
  post_density = density(resid_samples)
  post_density$y = post_density$y/length(resid_samples)
  post_mode = post_density$x[which.max(post_density$y)]
  xlim=max(post_density$x)
  colors=c('Prior'='blue','Posterior'='red')
  
  g <- ggplot(as.data.frame(model$VCV), aes(x=units)) + 
    geom_density(aes(y=..scaled..,color='Posterior'))+
    geom_line(aes(x = xvals, 
                  y = dinvgamma(xvals,shape = nu./2,scale = nu.*V_E/2),
                  color='Prior'))+
    geom_vline(aes(xintercept = prior_mode, color='Prior'))+
    geom_vline(aes(xintercept = post_mode, color='Posterior'))+
    labs(color='Legend')+
    ylab('density')+
    scale_color_manual(values=colors)+
    ggtitle('units')+
    annotate('text', x=c(prior_mode,post_mode),y=c(0.8*max(post_density$y),0.6*max(post_density$y)),label=c(as.character(round(prior_mode,4)),as.character(round(post_mode,4))),colour = c('blue','red'))
  print(g) 
}

#function to plot prior and posterior
plot_priorandpost_two <- function(model,prior,xlim=0){
  if(xlim==0) xlim=max(model$VCV)
  xvals = seq(0,xlim,length.out=nrow(model$VCV))  
  
  param_list = colnames(model$VCV)
  
  for(i in 1:length(param_list)){
    #get prior parameters for the current variance component
    if(i==length(param_list)) {prior_params = prior$R #last one is the residual
    } else prior_params = prior$G[[i]] #others are elements of G
    
    #generate prior density (y values to match the existing xvals vector)
    if(length(prior_params)==4){ #parameter-expanded priors
      prior_density = df(xvals/prior_params$alpha.V, df1 = 1, df2 = prior_params$nu, ncp = (prior_params$alpha.mu^2)/prior_params$alpha.V)
    } else if(length(prior_params)==2){ #regular inverse-gamma priors
      prior_density = dinvgamma(xvals,shape = prior_params$nu/2,scale = prior_params$nu*prior_params$V/2)
    } else print('Warning: cannot handle this type of prior')
    
    #get prior mode to use for plotting    
    prior_mode = xvals[which.max(prior_density)]
    if(prior_mode==xvals[length(xvals)]) print('Warning: xlim is too small')
    
    #grab posterior samples and convert to density
    post_samples = model$VCV[,i]
    post_density = density(post_samples)
    #post_density$y = post_density$y/post_density$n
    post_mode = post_density$x[which.max(post_density$y)]
    
    #do some scaling if necessary
    multfactor=1
    if((max(prior_density,na.rm=TRUE)*10) < max(post_density$y)){ #if the prior density is on a much smaller scale than the posterior
      multfactor = 0.4*(max(post_density$y)/max(prior_density,na.rm=TRUE))
      prior_density = prior_density * multfactor
    }
    
    colors=c('Prior'='blue','Posterior'='red')
    g <- ggplot(as.data.frame(model$VCV), aes(x=as.data.frame(model$VCV)[,i])) + 
      geom_density(aes(color='Posterior'))+
      xlim(0,xlim)+
      geom_line(aes(x = xvals, 
                    y = prior_density,
                    color='Prior'))+
      labs(color='Legend')+
      scale_color_manual(values=colors)+
      scale_y_continuous(name='Posterior density',sec.axis = sec_axis(~./multfactor, name='Prior density')) +
        theme(axis.text.y.right=element_text(color='blue'),
            axis.title.y.right=element_text(color='blue'),
            axis.title.y.left=element_text(color='red'),
            axis.text.y.left=element_text(color='red'),
            legend.position = 'none') + 
      xlab(paste(param_list[i],'variance'))+
      ggtitle(param_list[i])
    print(g)
    
  }
}
