
options(markdown.HTML.stylesheet = "path/to/a/custom/style.css")
library(DynamicCorrelation)
library(MASS)

set.seed(2019)
n = 30
k = 3
nullX = as.matrix(rnorm(30), ncol=1)
Sigma = matrix(0.2, 3, 3)
diag(Sigma) = 1.5
nullY = MASS::mvrnorm(n, rep(0,3), diag(3))

set.seed(2019)
n = 30
k = 3
Y = matrix(NA, n, k)
X = rnorm(30)
rho = X * 0.2
rho = rho - min(rho)
X = matrix(X, ncol=1)
for (i in 1:n){
  Sigma = matrix(rho[i], k, k)
  diag(Sigma) = 1
  Y[i,] = MASS::mvrnorm(1, rep(0,k), Sigma)
}

q = get_score(nullX, nullY[,1], nullY[,2])
p = 1-pchisq(q, 1)
print(paste("The score statistic is", q, "and the p-value is", p))

q = get_score(X, Y[,1], Y[,2])
p = 1-pchisq(q, 1)
print(paste("The score statistic is", q, "and the p-value is", p))

d = get_degree(as.matrix(X, ncol=1), Y[,1], Y[,2:3])
p = get_p_from_degree(Y[,1], Y[,2:3], d)
print(paste("The degree statistic is", d, "and the p-value is", p))

print(paste("The degree statistic is", d, "and the p-value is", p))



n = 30
k = 3
B = 1000
p = rep(NA, B)
for (b in 1:B){
  nullX = as.matrix(rnorm(30), ncol=1)
  Sigma = matrix(0.2, 3, 3)
  diag(Sigma) = 1.5
  nullY = MASS::mvrnorm(n, rep(0,3), Sigma)
  d = get_degree(as.matrix(nullX, ncol=1), nullY[,1], nullY[,2:3])
  p[b] = get_p_from_degree(Y[,1], Y[,2:3], d)
}


for (b in 1:B){
  nullX = as.matrix(rnorm(30), ncol=1)
  Sigma = matrix(0.2, 3, 3)
  diag(Sigma) = 1.5
  nullY = MASS::mvrnorm(n, rep(0,3), Sigma)
  q = get_newscores(nullX, nullY[,1], nullY[,2])
  p[b] = 1-pchisq(q, 1)
  # print(paste("The score statistic is", q, "and the p-value is", p))
}
print(sum(p<0.05))

