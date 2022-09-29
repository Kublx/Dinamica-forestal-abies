#########################################################
#################Modelo din√°mica poblacional Abies religiosa############
###crecimiento, mortalidad y regeneracion##############
############################################################
#files#
rodales_area <- read.csv(file="rodales_modelo_gobernanza.csv", header=T)
trees <- read.csv(file="numero_arboles_diamatreo_rodades.csv", header=T)
colnames(trees) <- c("Rodal","CD","Number")
########################
## Sampling functions ##
########################
subpop <- function(d = rodales_area) {
  res <-  d[sample(1:nrow(d), 1), ]
  numberrodal <- res[1]
  arearodal <- res[2]
  return(list(numberrodal=numberrodal,arearodal=arearodal))
  }


GetTrees <- function(T, Plots) {
  res <- NULL
  for (i in Plots) {
    res <- rbind(res, subset(T, Rodal %in% i))
  }
  return(res)
}

Sampled <- subpop ()
Plots <- Sampled$numberrodal

process<- function (Tree,rodales_area){
  Sampled <- subpop ()
  Plots <- Sampled$numberrodal
  Tree.List <- GetTrees (trees,Plots)
  stand <- Tree.List$Number 
  class <- Tree.List$CD
  basal_big_class <- 0.424455347
  BAB <- rep(0,100)
  TBA <- 3.142*(Tree.List[Tree.List[,2]>75,2]/200)^2###Partition the basal area of big trees >31 cm and add number of trees that the surplus of basal area represents
  BAB <- round(TBA/basal_big_class,digits=0)
  y <- sum(BAB)
  rodal <- c(0,stand[1:14])
  rodal[15] <- rodal[15]+y
  rodal[1]<- 30
  res <- rodal
  return (res)
}   
    

  




exe <- function(Y,Tree,rodales_area) {
  rodal<-process(Tree,rodales_area)
  n <- length(rodal)
  N1s <- rep(1, n)
  N0s <- rep(0, n)
  dbhl <- c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71) # Vector de limite inferior, clases diametricas
  dbhu <- c(6,11,16,21,26,31,36,41,46,51,56,61,66,71,76) # Vector de limite superior, clases diametricas
  dbhq <- sqrt((dbhu^3-dbhl^3)/((dbhu-dbhl) * 3)) # valor cuadr√É¬°tico de DBH para asignar un valor de rango a cada clase diam√É¬©trica
  baq  <- (dbhq/2)^2 * pi/1e4 # Area basal cuadratica 
  
  Rodal <- matrix(0, Y, n, byrow = TRUE)
  Recruits <- numeric(Y)
  BA <- numeric(Y)
  Senescence <- matrix(0, Y, n, byrow = TRUE)
  Crecimiento <- matrix(0, Y, n, byrow = TRUE)
  
  #para p. montezumae##knownpoints <- data.frame(x <- c(1,5.2,15.2,21,33.4,35,40.1,42,44.9),
                          ##  y <- c(0.08,1,2,1.8,1.1,0.8,0.5,0.2,0.2))
  
  #knownpoints <- data.frame(x <- c(1,3.5,10,17.3,30.2,41.5,50,70, 100),
                            #y <- c(1,1,1.4,1.38,1.33,1.1, 1.0,1.0,0.9))
  ##transition probabilities from Fir demography paper
  knownpoints <- data.frame(x <- c(1,3,5,10,20,25,30,35,40,45,50,60,100),
                            y <- c(0.0044,0.0491,0.0491,0.0294,0.0294,0.099,0.099,0.101,0.101,0.117,0.117,0.011,0.011))
  
  yos <- approx(knownpoints$x, knownpoints$y)##interpolacion
  diameter <- yos$x
  annualdiameter <- yos$y ##proveniente de las interopolaciones...constante
  #RegenLagTime <- 30   # assume it takes a 1-yr old natural seedling 30 years to grow to
  # a sapling with a DBH 2.0 cm (first diameter class) data from JosÈ L. Gallardo-Salazar1
  #, Dante A. RodrÌguez-Trejo1*, Salvador Castro-Zavala. 2019
  
  
  #iw <- dbhu-dbhl # clases de 5 cm
  ###Mortality TemiÒo-ViÒota et al. 2016. Modelling initial mortality of Abies religiosa
  ###produce altas mortalidades
  #a<- -(1.7033-0.0185*dbhq)
  #pmort <- 1/(1+exp(a))
  
  ###probabilidad de mortalidad de Abies (Vieilledent et al. 2009). 
 
  knownmortality <- data.frame(x <- c(2,5,15,45,75,100),
                            y <- c(0.05,0.043,0.0197,0.0029,0.0009,0.0006))
  yos_mort <- approx( knownmortality$x,  knownmortality$y)
  diameter_mort <- yos_mort$x
  mortality_prob <- yos_mort$y 
  #c1 <- 0.536#1-20 
  #c2 <- 0.157#20-40
  #c3 <- 0.371#>40 
  #probabilidad_mortalidad_1_20 <- runif(24000, c1-0.02,c1+0.02)
  #probabilidad_mortalidad_20_40 <- runif(24000, c2-0.02,c2+0.02)
  #probabilidad_mortalidad_40 <- runif(24000, c3-0.02,c3+0.02)
 

  
  
  for (y in 1:Y) {
    
    semillas <- floor(runif(2000,65,714)) #(#hacer ejercicio calibracion, 1500 a 2000 pl·ntulas por hectarea; Guzm·n-Aguilar et al. 2020)
    probabilityoftransitionseed_seedling <-0.0044
    plantulas <- round((sample(semillas,1,replace=TRUE))* probabilityoftransitionseed_seedling) 
    
    
    anndiam <- rep(0,n-1)
   
    diame <- rep(0,n-1)
    for (i in 1:length(anndiam)){
      vector <- which(diameter >= dbhq[i] & diameter<dbhq[i+1])
      randomvector <- sample(vector,1)
      tmp1 <- diameter[randomvector]
      tmp <- annualdiameter[randomvector] 
     
      #randomanndiam <- sample(runif(300, tmp-0.02,tmp+0.02),1)##aleatoriedad den el modelo
      anndiam[i] <- tmp
      diame[i] <- tmp1
    }
    
    vector15 <- which(diameter>=floor(dbhq[n-1]) & diameter<floor(dbhq[n]))
    ann15 <- sample(annualdiameter[vector15],1) 
    diam15 <- sample(diameter[vector15],1) 
    anndiam[15]<- ann15
    
    adi<- anndiam
    graduating <- adi
    #graduating <- 1/(iw/adi) # Probabilidad de graduacion
    rodal_crecimiento <- rep(0,length(rodal))
    for (x in 1:length(rodal_crecimiento)){
      if (rodal[x]==0){      
        rodal_crecimiento[x] <- 0
      }else{
        rodal_crecimiento[x] <- rbinom(1, rodal[x], graduating[x]) # Arboles que pasan a la siguiente clase
      }
    }
    
    
    rodal <- rodal - rodal_crecimiento # Resta de los arboles que se "graduan"
    rodal <- rodal + c(0, rodal_crecimiento[1:n-1]) # Suma de los arboles "graduados"
    #LastDCT <- rodal[length(rodal)]  # prevents loosing trees in last diameter class
    #rodal [length(rodal)] <- rodal[length(rodal)] + LastDCT
    BAIncremento <- sum(rodal_crecimiento * baq)    # patch BA increment due to growth
    # ##Reclutamiento
    #paper"Plantaciones forestales vs. regeneraci√≥n natural in situ: 
    #El caso de los pinos y la rehabilitaci√≥n en el Parque Nacional Cofre de Perote"
    ##Se estimaron de 1200 a 1600 pl√°ntulas por ha
    
    
    ##Mortalidad (size-dependent mortality rates)
    prob_mort <- rep(0,n-1)
    diame_mort <- rep(0,n-1)
    for (i in 1:length(prob_mort)){
      vector <- which(diameter_mort >= dbhq[i] &  diameter_mort<dbhq[i+1])
      tmp <- sample(mortality_prob[vector],1) 
      tmp1 <- sample( diameter_mort[vector],1) 
      #randomanndiam <- sample(runif(300, tmp-0.02,tmp+0.02),1)##aleatoriedad den el modelo
      prob_mort[i] <- tmp
      diame_mort[i] <- tmp1
    }
    
    vector_mort20 <- which(diameter_mort>=floor(dbhq[n-1]) & diameter_mort<floor(dbhq[n]))
    mort20 <- sample(mortality_prob[vector_mort20],1) 
    diam_mort20 <- sample(diameter_mort[vector_mort20],1) 
    prob_mort[n]<- mort20
    
    rodal_mortalidad <- rep(0,length(rodal))
    for (x in 1:length(rodal_mortalidad)){
      if (rodal[x]==0){      
        rodal_mortalidad[x] <- 0
      }else{
        rodal_mortalidad[x] <- rbinom(1, rodal[x], prob_mort[x]) # Arboles que mueren
      }
    }
    
    
    
    rodal <- rodal-rodal_mortalidad #solo por senescencia
    
    # Add recruitment
    #RegenCohorts <- rpois(RegenLagTime, sample(plantulas, 1, replace = TRUE))
    Recruits[y] <- plantulas
    #RegenCohorts <- c(NewRegen, RegenCohorts[1:RegenLagTime - 1])##pensar en esto cuando la regenarcion artificial
    rodal[1] <- rodal[1] + Recruits[y]
  
    ############ Salvar resultados simulacion##############
    Rodal[y, ] <- rodal # stand after regeneration, captures regeneration pulses
    Senescence[y, ] <- rodal_mortalidad
    Crecimiento[y, ] <-   rodal_crecimiento
  }
  
  res <- list(Rodal=Rodal, Senescence=Senescence,Crecimiento=Crecimiento, Recruits=Recruits)
  return(res)
}

set.seed(12)
exp <- exe(100)















#bsenescencia_1_20 <- rep(0, 10)
# for (a in 1:length(probsenescencia_1_20)){
#   rodaltmp1<- rodal[1:10]
#   if ( rodaltmp1[a]==0){
#     probsenescencia_1_20[a]<- 0
#   } else{
#     muertes_senescencia1 <- rbinom(1,  rodaltmp1[a], sample(probabilidad_mortalidad_1_20,1))
#     probsenescencia_1_20[a]<- muertes_senescencia1
#   }
# }
# 
# 
# probsenescencia_20_40 <- rep(0, 9)
# rodaltmp2<- rodal[11:19]
# for (b in 1:length(probsenescencia_20_40)){
#   if (rodaltmp2[b]==0){
#     probsenescencia_20_40[b]<- 0
#   }else{
#     muertes_senescencia2 <- rbinom(1,rodaltmp2[b], sample(probabilidad_mortalidad_20_40,1))
#     probsenescencia_20_40[b]<- muertes_senescencia2
#   }
# }
# 
# if (rodal[20]==0){      
#   probsenescencia_40 <- 0
# }else{
#   probsenescencia_40 <- rbinom(1, rodal[20], sample(probabilidad_mortalidad_40,1))
# }
# 
# senescencia <- c(probsenescencia_1_20,probsenescencia_20_40,probsenescencia_40 )
