
dataset  = read.spss("/home/anthony/Documentos/algor_genetico/tot2014s.sav", to.data.frame = TRUE)

clean_data = function(dados,beans){
  
  bean_size = round(length(dados)/beans)
  
  new_data = numeric(beans)
  j=0
  for(i in seq(1,length(dados),by=bean_size)){
    j=j+1
    new_data[j] = dados[i]
  }
  
  return(new_data)
}

datab = log(dataset$V2)

Y = hist(dataset$V2,100000,freq = TRUE)

func <- function(j,s,w){
  
  u = 1+(j/(s**2))
  
  p = ( ( (u-1)**u)/gamma(u) )*exp(-(u-1)/w)/(w**(1+u))
  
  dpdj = (p/(s**2))*( -log(w) - (2*log(abs(s))) + (u*(s**2)) + log(j) - (digamma(u)) - (1/(w)))
  
  dpds = (2*j*p/(s**3))*( (log(w)) + (2*log(abs(s))) - ((s**2)*u/j) + (digamma(u)) - (log(j)) + (1/j) ) 
  
  return(c(p,dpdj,dpds))
}

fitting <- function(j,s){  # x = [j,s,w]
  
  A = numeric(length(Y$density))
  B = numeric(length(Y$density))
  
  for(i in 1:length(Y$density)){
    

    funcVal = func(j,s,Y$mids[i])

    
    A[i] = -(Y$density[i] - funcVal[1])*funcVal[2]  # funcVal[1] = p(x), funcVal[2] =dp(x)/dj
    B[i] = -(Y$density[i] - funcVal[1])*funcVal[3]  # funcVal[3] =dp(x)/ds
    }
  return(sqrt((sum(A)**2) + (sum(B)**2)))
}

index = function(j,s){
  
  
  return(1 + (j/(s**2)))
  
}



GA = ga(type = "real-valued", 
        fitness =  function(x) -fitting(x[1],x[2]),
        lower = c(1e-5,1e-5), upper = c(5,5), 
        popSize = 50, maxiter = 1000, run = 100)

x1 <- x2 <- seq(0.001, 10, by = 0.1)
f <- outer(x1, x2, fitting)
persp3D(x1, x2, outer(x1,x2,f), theta = 50, phi = 20, color.palette = bl2gr.colors)

new_Y = numeric(length(Y$density))
for(i in 1:length(Y$density)){
  new_Y[i] = func(GA@solution[1],GA@solution[2],Y$mids[i])[1]
}

