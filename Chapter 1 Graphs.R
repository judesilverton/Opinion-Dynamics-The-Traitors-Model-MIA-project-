#########################################################

#Opinion dynamics model: Chapter 1

#########################################################


#SET WORKING DIRECTORY

setwd("~/Opinion dynamics")

#Define the function which we will use to calculate bounded confidence set I(i, x)

BCset <- function(x, i, epsilon, m){
  I <- numeric(m)
  for(j in 1:m){
    if(abs(x[i]-x[j]) <= epsilon){
      I[j] = 1
    }
  }
  return(I)
}



#########################################################
#Define the HK Model Function
########################################################

HKmodel <- function(m, maxtime, epsilon){
  
  
  #Generate opinions of m people between 0 and 1
  #Create empty list of m people 
  
  x <- runif(n=m, min=0, max=1)
  oplist <- vector("list", length = m)
  t <- 0
  
  
  #Put these random opinions in the list 
  
  for(i in 1:m){
    oplist[i] <- x[i]
  }
  
  #Perform the update rule until maxtime
  
  for(t in 1:maxtime){
    vec <- data.table::transpose(oplist)
    for(j in 1:m){
      y <- BCset(unlist(vec[t]), j, epsilon, m)
    oplist[[j]][t+1] <- (1/sum(y)) * sum(y*unlist(vec[t]))
    }
    
  }
  
  return(oplist)
  
}


#########################################################
#Run the model for epsilon = 0.30
#Corresponds to Figure 2 in Scholarly Report
########################################################

set.seed(20)

h1 <- HKmodel(40, 55, 0.3)


#Form the time vector for the x axis

t_vec <- 0

for(i in 1:55){
  t_vec[i+1] <- i
}


#Plot the results

pdf("Epsilon 30.pdf", width = 7, height = 6)

par(mar=c(3.7,4,2,2), mgp=c(2.4,0.7,0))
plot(t_vec, unlist(h1[1]), type = "l", col = "orange", lty = 1, lwd = 1.75, xlim = c(0, 15), ylim = c(0,1),
     xlab = "Time", ylab = "Opinion values", las = 1, cex.lab = 1.2)
title(expression(epsilon ~ "= 0.3, n = 40"))


#Plot all the opinion lines  over the time period
for(i in 2:40){
  lines(t_vec, unlist(h1[i]), type = "l", col ="orange", lwd=1.75)
}


dev.off()



#########################################################
#Run the model for epsilon = 0.15
#Corresponds to Figure 3 in Scholarly Report
########################################################

set.seed(38)

h2 <- HKmodel(40, 55, 0.15)

pdf("Epsilon 15.pdf", width = 7, height = 6)
par(mar=c(3.7,4,2,2), mgp=c(2.4,0.7,0))
v <- unlist(h2[1])

if(any(v >= 0.6) == TRUE){
  plot(t_vec, v, type = "l", col = "red", lty = 1, lwd = 1.7, xlim = c(0, 10), ylim = c(0,1),
       xlab = "Time", ylab = "Opinion values", las = 1, cex.lab = 1.2)
}

if(any(0.2 < v) == TRUE && any(v < 0.6) == TRUE){
  plot(t_vec, v, type = "l", col = "blue", lty = 1, lwd = 1.7, xlim = c(0, 10), ylim = c(0,1),
       xlab = "Time", ylab = "Opinion values", las = 1, cex.lab = 1.2)
}

if(any(0.2 >= v) == TRUE){
  plot(t_vec, v, type = "l", col = "orange", lty = 1, lwd = 1.7, xlim = c(0, 10), ylim = c(0,1),
       xlab = "Time", ylab = "Opinion values", las = 1, cex.lab = 1.2)
}

  
title(main = expression(epsilon ~ "= 0.15, n = 40"))

#Plot lines for each cluster

for(i in 2:40){
  v <- unlist(h2[i])
  if(any(v > 0.6) == TRUE){
  lines(t_vec, v, type = "l", col ="red", lwd = 1.7)
  }
  else if(any(0.2 < v) == TRUE && any(v < 0.6) == TRUE){
   lines(t_vec, v, type = "l", col ="blue", lwd = 1.7) 
  }
  else{
    lines(t_vec, v, type = "l", col ="orange", lwd = 1.7)
  }
}


dev.off()

