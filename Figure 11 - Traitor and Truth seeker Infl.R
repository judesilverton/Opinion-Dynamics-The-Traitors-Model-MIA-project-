#########################################################

#Opinion dynamics model: Chapter 4, Figure 11

########################################################

#THIS MODEL INCLUDES INFLUENCE LEVELS


library(DescTools)
library(raster)
library(ggplot2)
library(rasterVis)

#SET WORKING DIRECTORY

setwd("~/Opinion dynamics")


#Update rule WITH INFLUENCE LEVELS for opinion follower 

update <- function(x, i, epsilon, m, t, r){
  I <- numeric(m*m)
  SUM <- 0
  r[i] <- 1
  for(j in 1:m){
    #For each agent j, check if agent i and j's opinions are within epsilon, wrt standard norm
    if(norm(x[[i]][(m*t+1):(m*(t+1))]-x[[j]][(m*t+1):(m*(t+1))], type="2") <= epsilon[i]){
      I[((j-1)*m +1):(j*m)] <- rep(r[j],m)
      SUM <- SUM + (x[[j]][(m*t+1):(m*(t+1))] * I[((j-1)*m +1):(j*m)])
    }
    new <- SUM/(sum(I)/m)
  }
  return(new)
}

#Update rule WITH INFLUENCE LEVELS for truth seeker

updateTRUTHSEEKER <- function(x, i, epsilon, m, t, alpha, truth, r){
  I <- numeric(m*m)
  SUM <- 0
  r[i] <- 1
  for(j in 1:m){
    if(norm(x[[i]][(m*t+1):(m*(t+1))]-x[[j]][(m*t+1):(m*(t+1))], type="2") <= epsilon[i]){
      I[((j-1)*m +1):(j*m)] <- rep(r[j],m)
      SUM <- SUM + (x[[j]][(m*t+1):(m*(t+1))] * I[((j-1)*m +1):(j*m)])
    }
    new <- alpha*truth + (1-alpha)*SUM/(sum(I)/m)
  }
  return(new)
}


###########################################################
#Define our extended HK model function WITH INFLUENCE
#########################################################


HKmodel <- function(m, maxtime, epsilon, k, alpha, truth, lie, r){
  
  
  #generate opinions of m people 
  #create empty list of m people 
  
  x <- vector("list", length = m)
  
  #agent 1 is traitor
  #agent 2 is a truth seeker
  
  x[[1]][1:m] <- lie
  

  #Make a vector of probabilities (adding to 1), for agents 2 to m
  for(j in 2:m){
    #generate 5 random probabilities 
    a <- runif(n=k, min=0, max=1)
    P <- numeric(m)
    #generate k random agents for those probabilities to be of
    random_agents <- sample.int(m,k)
    P[random_agents] <- a
    P[j] <- 0
    P2 <- P/sum(P)
    x[[j]][1:m] <- P2
  }
  
  diff <- 0
  final <- 0
  
  #perform the update rule until maxtime
  
  for(t in 1:(maxtime-1)){ 
    x[[1]][(m*t+1):(m*(t+1))] <- update(x, 1, epsilon, m, t-1, r)
    x[[2]][(m*t+1):(m*(t+1))] <- updateTRUTHSEEKER(x, 2, epsilon, m, t-1, alpha, truth, r)
    diff <- diff + norm(x[[2]][(m*t+1):(m*(t+1))] - x[[2]][(m*(t-1)+1):(m*t)], type="2")
    count <- t
    for(i in 3:m){
      x[[i]][(m*t+1):(m*(t+1))] <- update(x, i, epsilon, m, t-1, r)
      diff <- diff + norm(x[[i]][(m*t+1):(m*(t+1))] - x[[i]][(m*(t-1)+1):(m*t)], type="2")
      
    }
    
    #Stop simulation if something close to a steady state has been reached (see report)
    
    if(diff < 10^(-3)){
      for(k in 1:m){
        max <- which.max(x[[k]][(m*t+1):(m*(t+1))][-k])
        if(max >= k){
          final[k] <- max+1 
        }
        else{
          final[k] <- max
        }
      }
      
      break
     
    }
    else{
      diff <- 0
    }
  }
  
  if(count == (maxtime-1)){
    for(k in 1:m){
      max <- which.max(x[[k]][(m*(maxtime-1)+1):(m*maxtime)][-k])
      if(max >= k){
        final[k] <- max+1 
      }
      else{
        final[k] <- max
      }
      
    } 
    } 

  
  return(final)
  
}


#Set up inputs for simulations (truth value, traitor strategy)

#Truth value of agent 1 being a traitor

traitor <- numeric(10)
traitor[1] <- 1

#Opinion "fully suspecting" agent 2 (truth seeker)

lie <- numeric(10)
lie[2] <- 1


#########################################################
#Run simulations
########################################################

set.seed(2208)

iters <- 500

sims <- numeric(iters)
max <- numeric(15*20)

#Influence level for traitor will be varied, all others set to 1

infl <- rep(1,10)

START <- Sys.time()
for(n in 1:15){
  for(u in 1:20){
  eps_n <- rep((1/10) * n  , 10)
  eps_n[1] <- 0
  infl[1] <- u
  
for(l in 1:iters){
  res <- HKmodel(10, 50, eps_n, 4, 1/2, traitor, lie, infl)
  votes <- Mode(res)
  out <- votes[sample.int(length(votes), size = 1)]
  sims[l] <- out
    
}
  height <- table(unlist(sims))
  max[20*(n-1) + u] <- which.max(height)
  
  }
  print(n)
}
END <- Sys.time()

END - START

#Store output

write.csv(max, "~/Opinion dynamics/Figure 11.csv")

maxxx <- max[1:300]

#Rearrange so that graph is plotted correctly

maxxx <- maxxx[c(281:300, 261:280, 241:260, 221:240, 201:220, 181:200, 161:180, 141:160, 121:140, 101:120, 81:100, 61:80, 41:60, 21:40, 1:20)]

#Sort agents into categories

maxxx[maxxx>2] <- 3

#Sort data into format ready to plot

r <- raster(xmn=0, xmx=20, ymn =0, ymx=1.5, nrows = 15, ncols = 20)
r[] <- maxxx
r_points <- rasterToPoints(r)
r_df <- data.frame(r_points)

#Plot Figure 11

pdf("Figure 11.pdf", width =9, height = 6)
ggplot(r_df, aes(x, y, fill= as.factor(layer))) + 
  geom_tile(color = "black") + scale_fill_manual(breaks=c(0,1,2,3), labels=c("","Traitor","Truth seeker","Other"),
                                                 values=c("white","darkred","slateblue3","palegreen1")) + labs(x = "Influence level of traitor", y = expression(epsilon)) +
  theme(axis.title.y = element_text(angle=0, vjust=0.46, margin = margin(r=20)), axis.text = element_text(size=14),
        axis.title = element_text(size=15), axis.title.x = element_text(angle=0, margin = margin(t=16)),
        legend.key.size = unit(1, 'cm'), legend.title = element_text(size=14), 
        legend.text = element_text(size=13)) + labs(fill='Banished contestant')+
  scale_y_continuous(breaks=seq(0,1.5, by = 0.1)) + scale_x_continuous(breaks=seq(0,20, by = 2))
dev.off()



