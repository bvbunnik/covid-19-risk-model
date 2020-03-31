rm(list=ls())

library("deSolve"); library("ggplot2"); library("cowplot"); library(tidyverse)

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}


SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSp = -beta_pp*Sp*Ip - beta_ps*Sp*Is
    dSs = -beta_ss*Ss*Is - beta_sp*Ss*Ip
    dIp = beta_pp*Sp*Ip + beta_sp*Sp*Is - gamma*Ip
    dIs = beta_ss*Ss*Is + beta_ps*Ss*Ip - gamma*Is
    dR = gamma*Ip + gamma*Is
    return(list(c(dSp, dSs, dIp,dIs,dR)))
  })
}



init <- c(Sp = (0.85-0.001),Ss=0.15, Ip = 0.001, Is=0, R = 0)
times <- seq(0,730,by = 1)
R0=2.4
gamma = 1/(GenTime(4.6,2.4))

parms = c(gamma,
          beta_pp = gamma*R0,
          beta_ps = gamma*R0,
          beta_ss = gamma*R0,
          beta_sp = gamma*R0)
out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))

out2 <- out1 %>% pivot_longer(-time, names_to = "compartment", values_to = "fraction")

ggplot(out2, aes(x=time, y=fraction, colour=compartment)) + geom_line(size=1.05)
