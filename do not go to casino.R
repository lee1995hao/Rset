resu <- c()
for (i in 1:500000) {
  resu[i]<- sample(1:6,1)
}
r <- table(resu)
abs(r[1] - r[2])/5000000
mean(r)
hist(resu)