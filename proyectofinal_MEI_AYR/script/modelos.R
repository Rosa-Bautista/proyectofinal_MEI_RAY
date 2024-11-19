#### 


library(deSolve)
library(ggplot2)
library(gridExtra)


## MODELO 1 PROPUESTO

sin.tratamiento <- function(time,state,parameters){
  with(as.list(c(state, parameters)), {
    dSm <- -beta1 * Sm * Ih - mu1 * Sm + nu1 *(Sm + Em + Im)
    dEm <- beta1 * Sm * Ih - sigma * Em - mu1 * Em
    dIm <- sigma * Em - mu1 * Im
    dSh <- nu2 * (Sh + Eh + Ih) - beta2 * Sh * Im - mu2 * Sh
    dEh <- beta2 * Sh * Im - alpha * Eh - mu2 * Eh
    dIh <- alpha * Eh - mu2 * Ih
    return(list(c(dSm, dEm, dIm, dSh, dEh, dIh)))
  })
}


##parametros
b1 <- 0.048
b2 <- 24.96
nu1 <- 0.071
nu2 <- 0.004
mu1 <- 0.071
mu2 <- 0.004
s <- 0.0714
a <- 0.00236

parameters <- c(beta1= b1, beta2= b2, nu1= nu1, nu2= nu2,
                mu1= mu1, mu2= mu2, sigma= s, alpha= a)
state <- c(Sm= 47, Em= 0, Im= 53, Sh= 72, Eh= 0, Ih= 253)
time <- seq(0,1000, by= 1)

out_modelo1 <- ode(y = state, times = time, func = sin.tratamiento, parms = parameters)

##graficas

matplot(out_modelo1 [,1], out_modelo1 [, 5:7], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Susceptible", "Expuestos","Infectados"), col = 1:3, lty=1)

matplot(out_modelo1 [,1], out_modelo1 [, 2:4], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Suscepetible", "Expuestos","Infectados"), col = 1:3, lty=1)

## MODELO 2 PROPUESTO
tratamiento.infectivos <- function(time,state,parameters){
  with(as.list(c(state, parameters)), {
    dSm <- -beta1 * Sm * Ih - mu1 * Sm + nu1 *(Sm + Em + Im) 
    dEm <- beta1 * Sm * Ih - sigma * Em - mu1 * Em
    dIm <- sigma * Em - mu1 * Im
    dSh <- nu2 * (Sh + Eh + Ih + Th) - beta2 * Sh * Im - mu2 * Sh + rho * Th
    dEh <- beta2 * Sh * Im - alpha * Eh - mu2 * Eh
    dIh <- alpha * Eh - gamma * Ih - mu2 * Ih
    dTh <- gamma * Ih - rho * Th - mu2 * Th
    return(list(c(dSm, dEm, dIm, dSh, dEh, dIh, dTh)))
  })
}

##parametros
b1 <- 0.048
b2 <- 24.96
nu1 <- 0.071
nu2 <- 0.004
mu1 <- 0.071
mu2 <- 0.004
s <- 0.0714
a <- 0.00236
g <- 0.2
r <- 0.00018

parameters <- c(beta1= b1, beta2= b2, nu1= nu1, nu2= nu2,
                mu1= mu1, mu2= mu2, sigma= s, alpha= a,
                gamma= g, rho= r)
state <- c(Sm= 47, Em= 0, Im= 53, Sh= 72, Eh= 0, Ih= 253, Th= 0)
time <- seq(0,5478, by= 1)

out_modelo2 <- ode(y = state, times = time, func = tratamiento.infectivos, parms = parameters)

##graficas

matplot(out_modelo2 [,1], out_modelo2 [, 5:8], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Susceptible", "Expuestos","Infectados", "Tratamiento"), col = 1:4, lty=1)

matplot(out_modelo2 [,1], out_modelo2 [, 2:4], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Suscepetible", "Expuestos","Infectados"), col = 1:3, lty=1)

## MODELO 3 PROPUESTO
tratamiento.aleatorio <- function(time,state,parameters){
  with(as.list(c(state, parameters)), {
    dSm <- -beta1 * Sm * Ih - mu1 * Sm + nu1 *(Sm + Em + Im) 
    dEm <- beta1 * Sm * Ih - sigma * Em - mu1 * Em
    dIm <- sigma * Em - mu1 * Im
    dSh <- nu2 * (Sh + Eh + Ih + Th + Mh) - beta2 * Sh * Im - mu2 * Sh + rho * Th - gamma * Sh * pi + theta * Mh
    dEh <- beta2 * Sh * Im - alpha * Eh - mu2 * Eh - gamma * Eh * pi
    dIh <- alpha * Eh - gamma * pi * Ih - mu2 * Ih
    dTh <- gamma * pi * Ih - rho * Th - mu2 * Th
    dMh <- gamma * pi * (Sh + Eh) - theta * Mh - mu2 * Mh
    return(list(c(dSm, dEm, dIm, dSh, dEh, dIh, dTh, dMh)))
  })
}


##parametros
b1 <- 0.048
b2 <- 24.96
nu1 <- 0.071
nu2 <- 0.004
mu1 <- 0.071
mu2 <- 0.004
s <- 0.0714
a <- 0.00236
g <- 0.2
r <- 0.00018
p <- 0.8
t <- 0.0055

parameters <- c(beta1= b1, beta2= b2, nu1= nu1, nu2= nu2,
                mu1= mu1, mu2= mu2, sigma= s, alpha= a,
                gamma= g, rho= r, pi= p, theta= t)
state <- c(Sm= 47, Em= 0, Im= 53, Sh= 72, Eh= 0, Ih= 253, Th= 0, Mh= 0)
time <- seq(0,5478, by= 1)

out_modelo3 <- ode(y = state, times = time, func = tratamiento.aleatorio, parms = parameters)

##graficas

matplot(out_modelo3 [,1], out_modelo3 [, 5:9], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Susceptible", "Expuestos","Infectados", "Tratamiento", "Inmunes temp"), col = 1:5, lty=1)

matplot(out_modelo3 [,1], out_modelo3 [, 2:4], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Suscepetible", "Expuestos","Infectados"), col = 1:3, lty=1)

## MODELO 4 PROPUESTO
tratamiento.general <- function(time,state,parameters){
  with(as.list(c(state, parameters)), {
    dSm <- -beta1 * Sm * Ih - mu1 * Sm + nu1 *(Sm + Em + Im) 
    dEm <- beta1 * Sm * Ih - sigma * Em - mu1 * Em
    dIm <- sigma * Em - mu1 * Im
    dSh <- nu2 * (Sh + Eh + Ih + Th + Mh) - beta2 * Sh * Im * (1 - (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh)))) - mu2 * Sh + rho * Th - gamma * Sh * (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh))) + theta * Mh
    dEh <- beta2 * Sh * Im * (1 - (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh)))) - alpha * Eh * (1 - (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh)))) - mu2 * Eh - gamma * Eh * (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh)))
    dIh <- alpha * Eh * (1 - (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh)))) - gamma * Ih - mu2 * Ih
    dTh <- gamma * Ih - rho * Th - mu2 * Th
    dMh <- gamma * (pi - ((Ih + Th)/(Sh + Eh + Ih + Th + Mh))) * (Sh + Eh) - theta * Mh - mu2 * Mh
    return(list(c(dSm, dEm, dIm, dSh, dEh, dIh, dTh, dMh)))
  })
}


##parametros
b1 <- 0.00048
b2 <- 24.96
nu1 <- 0.071
nu2 <- 0.004
mu1 <- 0.071
mu2 <- 0.004
s <- 0.0714
a <- 0.00236
g <- 0.2
r <- 0.00018
p <- 0.85
t <- 0.0055

parameters <- c(beta1= b1, beta2= b2, nu1= nu1, nu2= nu2,
                mu1= mu1, mu2= mu2, sigma= s, alpha= a,
                gamma= g, rho= r, pi= p, theta= t)
state <- c(Sm= 47, Em= 0, Im= 53, Sh= 72, Eh= 0, Ih= 253, Th= 0, Mh= 0)
time <- seq(0,5478, by= 1)

out_modelo4 <- ode(y = state, times = time, func = tratamiento.general, parms = parameters)

##graficas

matplot(out_modelo4 [,1], out_modelo4 [, 5:9], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Susceptible", "Expuestos","Infectados", "Tratamiento", "Inmunes temp"), col = 1:5, lty=1)

matplot(out_modelo4 [,1], out_modelo4 [, 2:4], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Suscepetible", "Expuestos","Infectados"), col = 1:3, lty=1)

matplot(out_modelo4 [2000:5000,1], out_modelo4 [2000:5000, c(5, 7)], type = "l", xlab= "tiempo", ylab="poblacion", lwd= 2, lty= 1)
legend("topright", c("Suscepetible", "Infectados"), col = 1:2, lty=1)

