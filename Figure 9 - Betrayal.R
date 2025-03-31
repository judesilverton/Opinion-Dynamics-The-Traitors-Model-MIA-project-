#########################################################

#Opinion dynamics model: Chapter 4, Figure 9 BETRAYAL

########################################################

library(DescTools)
library(raster)
library(ggplot2)
library(rasterVis)

#SET WORKING DIRECTORY

setwd("~/Opinion dynamics")


#Update rule for opinion averager

update <- function(x, i, epsilon, m, t){
  I <- numeric(m*m)
  SUM <- 0
  for(j in 1:m){
    #For each agent j, check if agent i and j's opinions are within epsilon, wrt standard norm
    if(norm(x[[i]][(m*t+1):(m*(t+1))]-x[[j]][(m*t+1):(m*(t+1))], type="2") <= epsilon[i]){
      I[((j-1)*m +1):(j*m)] <- rep(1,m)
      SUM <- SUM + (x[[j]][(m*t+1):(m*(t+1))] * I[((j-1)*m +1):(j*m)])
    }
    new <- SUM/(sum(I)/m)
  }
  return(new)
}

#Update rule for truth seeker

updateTRUTHSEEKER <- function(x, i, epsilon, m, t, alpha, truth){
  I <- numeric(m*m)
  SUM <- 0
  for(j in 1:m){
    if(norm(x[[i]][(m*t+1):(m*(t+1))]-x[[j]][(m*t+1):(m*(t+1))], type="2") <= epsilon[i]){
      I[((j-1)*m +1):(j*m)] <- rep(1,m)
      SUM <- SUM + (x[[j]][(m*t+1):(m*(t+1))] * I[((j-1)*m +1):(j*m)])
    }
    new <- alpha*truth + (1-alpha)*SUM/(sum(I)/m)
  }
  return(new)
}






###########################################################
#Define our extended HK model function
#########################################################

HKmodel <- function(m, maxtime, epsilon, k, alpha, truth, truth2, lie){
  
  
  #generate opinions of m people 
  #create empty list of m people 
  
  x <- vector("list", length = m)
  
  #agent 1,2 are traitors
  #agent 3,4 are truth seekers
  
  x[[1]][1:m] <- lie
  
  

  #Make a vector of probabilities (adding to 1), for agents 2 to m
  for(j in 2:m){
    #generate k random probabilities 
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
    x[[1]][(m*t+1):(m*(t+1))] <- update(x, 1, epsilon, m, t-1)
    x[[2]][(m*t+1):(m*(t+1))] <- updateTRUTHSEEKER(x, 2, epsilon, m, t-1, 1/4, truth)
    x[[3]][(m*t+1):(m*(t+1))] <- updateTRUTHSEEKER(x, 3, epsilon, m, t-1, alpha, truth)
    x[[4]][(m*t+1):(m*(t+1))] <- updateTRUTHSEEKER(x, 4, epsilon, m, t-1, alpha, truth2)
    diff <- diff + norm(x[[2]][(m*t+1):(m*(t+1))] - x[[2]][(m*(t-1)+1):(m*t)], type="2")
    diff <- diff + norm(x[[3]][(m*t+1):(m*(t+1))] - x[[3]][(m*(t-1)+1):(m*t)], type="2")
    diff <- diff + norm(x[[4]][(m*t+1):(m*(t+1))] - x[[4]][(m*(t-1)+1):(m*t)], type="2")
    count <- t
    for(i in 5:m){
      x[[i]][(m*t+1):(m*(t+1))] <- update(x, i, epsilon, m, t-1)
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

#Truth value of agent 2 being a traitor

traitor2 <- numeric(10)
traitor2[2] <- 1

#Opinion "fully suspecting" agent 3

lie <- numeric(10)
lie[3] <- 1



#########################################################
#Run simulations
########################################################

set.seed(2048)

iters <- 500

sims <- numeric(iters)
max <- numeric(10*10)

START <- Sys.time()
for(n in 1:10){
  for(u in 1:10){
  eps_n <- rep((1/20) * n, 10)
  eps_n[1] <- 0
  alpha_u <- (1/20)*u 
for(l in 1:iters){
  res <- HKmodel(10, 50, eps_n, 4, alpha_u, traitor, traitor2, lie)
  votes <- Mode(res)
  out <- votes[sample.int(length(votes), size = 1)]
  sims[l] <- out
    
}
  height <- table(unlist(sims))
  max[10*(n-1) + u] <- which.max(height)
  
  }
  print(n)
}
END <- Sys.time()

END - START

#Store output

write.csv(max, "~/Opinion dynamics/Figure 9 BETRAYAL.csv")


maxxx <- max[1:100]


#Rearrange so that graph is plotted correctly


maxxx <- maxxx[c(91:100, 81:90, 71:80, 61:70, 51:60, 41:50, 31:40, 21:30, 11:20, 1:10)]

maxxx[maxxx>4] <- 5

#Sort data into format ready to plot

r <- raster(xmn=0, xmx=0.5, ymn =0, ymx=0.5, nrows = 10, ncols = 10)
r[] <- maxxx
r_points <- rasterToPoints(r)
r_df <- data.frame(r_points)


#Plot Figure 9

pdf("Figure 9 BETRAYAL.pdf", width =8, height = 6)
ggplot(r_df, aes(x, y, fill= as.factor(layer))) + 
  geom_tile(color = "black") + scale_fill_manual(breaks=c(0,1,2,3,4, 5), labels=c("","Traitor 1", "Traitor 2", "Truth seeker 1", "Truth seeker 2", "Other"),
                                                 values=c("white","darkred","indianred3", "slateblue3","cadetblue3", "palegreen1")) + labs(x = expression(alpha), y = expression(epsilon)) +
  theme(axis.title.y = element_text(angle=0, vjust=0.46, margin = margin(r=20)), axis.text = element_text(size=14),
        axis.title = element_text(size=19), axis.title.x = element_text(angle=0, margin = margin(t=16)),
        legend.key.size = unit(1, 'cm'), legend.title = element_text(size=14), 
        legend.text = element_text(size=13)) + labs(fill='Banished contestant')+
  scale_y_continuous(breaks=seq(0,0.5, by = 0.1)) + scale_x_continuous(breaks=seq(0,0.5, by = 0.1))
dev.off()

