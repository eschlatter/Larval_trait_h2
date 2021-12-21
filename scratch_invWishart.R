library(extraDistr) #has inverse gamma
library(MCMCpack) #for inverse gamma and inverse wishart distributions

#prior
prior_invgamma_params <- list(V = 1,nu = 0.000001)
param_list <- prior_invgamma_params
nu <- param_list$nu
V <- param_list$V
#corresponding shape and scale (invgamma) parameters
shape=nu/2
scale=nu*V/2

xrange <- seq(0.0001,1,.0001)

#calculate inverse gamma
y_invgamma=MCMCpack::dinvgamma(xrange,shape,scale)
ymax=max(y_invgamma)

#calculate inverse wishart
distys <- data.frame(x=xrange,y=NA,dist='invgamma')
for(i in 1:nrow(distys)){distys$y[i]=diwish(distys$x[i],V,nu)}
distys_scaled <- ymax*distys$y/max(distys$y) #scale it to have the same maximum as the inverse gamma version


ggplot()+
  #geom_line(aes(x=xrange, y=extraDistr::dinvgamma(xrange,alpha=shape,beta=scale)))+ #use extraDistr invgamma function
  geom_line(aes(x=xrange,y=y_invgamma))#+ #equivalently, MCMCpack invgamma function
#  geom_line(aes(x=xrange,y=distys_scaled),color='blue') #using diwish from MCMCpack -- NOT equivalent. much steeper.
