library(DynamicCorrelation)
source("R/utilities.R")
sourceCpp("src/utilities.cpp")
set.seed(2019)
n = 100
k = 3


# When there is nonlinearity unaccounted for
X1 = as.matrix(scale(rnorm(n))*sqrt(n)/sqrt(n-1), ncol=1)
X2 = as.matrix(scale(rnorm(n))*sqrt(n)/sqrt(n-1), ncol=1)

q = p = rep(NA, 1000)
for (sim in 1:1000){
  nullnoise = MASS::mvrnorm(n, rep(0,3), diag(3))

  Y1 = X1^2 + nullnoise[,1]
  Y2 = X1^2 + nullnoise[,2]
  Y3 = X1^2 + nullnoise[,3]

  r1 = scale(resid(lm(Y1~X1)))*sqrt(n)/sqrt(n-1)
  r2 = scale(resid(lm(Y2~X1)))*sqrt(n)/sqrt(n-1)
  # r3 = scale(resid(lm(Y1~X1)))*sqrt(n)/sqrt(n-1)

  q[sim]= get_score(X1, r1, r2)
  p[sim] = 1-pchisq(q[sim], 1)
  # print(paste("The score statistic is", q, "and the p-value is", p))
}
print(sum(p<0.05)) ;hist(p)



# When there's a confounder that is correlated with both X and Y
X1 = as.matrix(scale(rnorm(n))*sqrt(n)/sqrt(n-1), ncol=1)
X2 = as.matrix(scale(X1+rnorm(n))*sqrt(n)/sqrt(n-1), ncol=1)
q2 = p2 = rep(NA, 1000)
for (sim in 1:1000){
  nullnoise = MASS::mvrnorm(n, rep(0,3), diag(3))

  Y1 = X1 + X2 + nullnoise[,1]
  Y2 = X1 + X2 + nullnoise[,2]
  Y3 = X1 + X2 + nullnoise[,3]

  r1 = scale(resid(lm(Y1~X1)))*sqrt(n)/sqrt(n-1)
  r2 = scale(resid(lm(Y2~X1)))*sqrt(n)/sqrt(n-1)
  # r3 = scale(resid(lm(Y1~X1)))*sqrt(n)/sqrt(n-1)

  q2[sim]= get_score_c(X1, r1, r2)
  p2[sim] = 1-pchisq(q2[sim], 1)
  # print(paste("The score statistic is", q, "and the p-value is", p))
}
fdr2 = sum(p2<0.05)
print(fdr2); hist(p2)
