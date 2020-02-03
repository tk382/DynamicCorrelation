Rcpp::sourceCpp("src/utilities.cpp")

get_newgrad = function(x, y1, y2){
  a = mean(y1*y1)
  b = mean(y2*y2)
  d = mean(y1*y2)
  return(sum((y1*y2-d)*x))
}

get_newfisher = function(x, y1, y2){
  a = mean(y1*y1)
  b = mean(y2*y2)
  d = mean(y1*y2)
  xmat = cbind(rep(1, length(x)),x)
  tmpmat = solve(t(xmat) %*% xmat)
  return(tmpmat)
}

get_newscores = function(x, y1, y2, option = 1){
  a = mean(y1*y1)
  b = mean(y2*y2)
  d = mean(y1*y2)
  const = 1
  const = (a*b + d^2) / (a*b-d^2) / (a*b-d^2)
  const2 = (a*b + d^2)^3 / (a*b-d^2)^4
  const3 = 1/(a*b+d^2)
  if(option == 1){
    return(const *
             get_newgrad(x,y1,y2)^2 /length(x))
  }else if(option == 2){
    return(const2 *
             get_newgrad(x,y1,y2)^2 /length(x))
  }else if (option == 3){
    return(const3 *
             get_newgrad(x, y1, y2)^2 /length(x))
  }

}

#### null scenario ####
# set.seed(2019)
n = 100
k = 2

B = 1000
newscores = newscores2 = newscores3 = rep(NA, B)
scores = rep(NA, B)
abd = matrix(NA, B, 4)
for (j in 1:B){
  nullX = as.matrix(rnorm(n), ncol=1)
  nullX = scale(nullX) * sqrt(n) / sqrt(n-1)
  nullY = MASS::mvrnorm(n, rep(0,2), diag(2))
  nullY = scale(nullY) * sqrt(n) / sqrt(n-1)
  newscores[j] = get_newscores(nullX, nullY[,1], nullY[,2])
  newscores2[j] = get_newscores(nullX, nullY[,1], nullY[,2], option = 2)
  newscores3[j] = get_newscores(nullX, nullY[,1], nullY[,2], option = 3)
  scores[j] = get_score(nullX, nullY[,1], nullY[,2])
  abd[j,] = c(sum(nullY[,1]^2),
              sum(nullY[,2]^2),
              sum(nullY[,1]*nullY[,2]),
              sum(nullY[,1]*nullY[,2]*nullX))
}
par(mfrow = c(2,2))
ref = qchisq(seq(0,1,length=1000), 1)
plot(density(scores), ylim = c(0,1), xlim = c(0, 4), col = 'green')
lines(density(newscores3), col = 'red')
lines(density(ref), col = 'black')
plot(newscores3 ~ scores); abline(0,1,col = 'red')
qqplot(ref, newscores3); abline(0,1,col = 'red')
qqplot(ref, scores); abline(0,1,col = 'red')
print(c(mean(newscores3), mean(scores)))
print(c(sum(newscores > qchisq(0.95, 1)), sum(scores>qchisq(0.95, 1))))
# hist(pchisq(newscores, 1)); hist(pchisq(scores, 1))


par(mfrow = c(1,3))

plot(newscores ~ abs(abd[,4]))
plot(newscores2 ~ abs(abd[,4]))
plot(newscores3 ~ abs(abd[,4]))

#
#### alternative ####
# set.seed(2019)
n = 100
k = 3

Anewscores = Ascores = rep(0, B)
for (j in 1:B){
  Y = matrix(NA, n, k)
  X = rnorm(n)
  rho = X * 0.1
  rho = rho - min(rho)
  X = matrix(X, ncol=1)
  for (i in 1:n){
    Sigma = matrix(rho[i], k, k)
    diag(Sigma) = 1
    Y[i,] = MASS::mvrnorm(1, rep(0,k), Sigma)
  }
  Anewscores[j] = get_newscores(X, Y[,1], Y[,2], option = 3)
  Ascores[j] = get_score(X, Y[,1], Y[,2])
}
plot(density(Ascores), ylim = c(0,1.5))
lines(density(Anewscores), col = 'red')
sum(Anewscores > qchisq(0.95, 1))
sum(Ascores > qchisq(0.95, 1))
