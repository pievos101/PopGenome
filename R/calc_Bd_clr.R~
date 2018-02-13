calc_Bd_clr <- function(D, global.mat, lambda){



D             <- D[!is.na(D)]
#D[is.na(D)] <- 0


n.valid.sites <- length(D)

if(length(D)==0){return(list(NaN,NaN))}
if(length(D)==1){return(list(NaN,NaN))}


#D         <- abs(D)
#D         <- round(D,digits=3)


# local Like
local.mat <- calc_Bd_CLR_table(D, lambda)


SAVECLR      <- numeric(100)

for  ( yy in 1:100){

#perm <- sample(1:(length(D)-1),length(D)-1,replace=FALSE)

perm <- sample(1:(length(D)-1),5,replace=TRUE)
perm <- sort(perm)


D.names      <- as.numeric(rownames(local.mat))
val          <- rep(NaN,length(D)-1) 
pos 	     <- rep(FALSE,length(D)-1)
neg	     <- rep(FALSE,length(D)-1)	

#Bd_local_vec <- rep(NaN,length(D)-1) 
#perm  <- 1:(length(D)-1)

for  ( xx in perm){

	row      <- match(D[xx],D.names)
	col      <- match(D[xx+1],D.names)


	if((D[xx]+D[xx+1])<=0){
		neg[xx] <- TRUE
	}

	if((D[xx]+D[xx+1])>=0){
		pos[xx] <- TRUE
        }

	#xxx      <- abs(D[xx]-D[xx+1])
	#row      <- match(xxx,D.names)

	#print("D1")
	#print(D[xx])

	#print("D2")
	#print(D[xx+1])


	#print("row")
	#print(row)
	#print("col")
	#print(col)
	
	#valx    <-  local.mat[row,col]/sum(local.mat)
	#valx    <-  sum(local.mat[row,])/sum(local.mat)
	#valxx   <-  sum(local.mat[,col])/sum(local.mat)
	#valy    <-  local.mat[row,col]/(sum(local.mat[row,])+sum(local.mat[,col]))	
       
	valz     <-  local.mat[row,col]/sum(local.mat) 
		
	#valz    <-  valz * (2-abs(mean(D,na.rm=TRUE)))

	#val[xx]  <- (valx + valxx) * valy * abs(D[xx]+D[xx+1]) #(1+abs(D[xx]+D[xx+1])) 	
	val[xx]  <-  valz #* (1-(abs(D[xx]+D[xx+1])/2)) #* (((abs(mean(D)))))  #* (1-abs(D[xx])) #+D[xx+1])) #*(2-abs(D[xx]-D[xx+1])) )	

	#print("value1")
	#print(value1)
	#print("value2")
	#print(value2)
	#print("val")
	#print(val[xx])

}



#DD      <- 0.0001 + (1-abs(D))
#Bd_local_vec[yy] <- sum(log(val[val>0]),na.rm=TRUE) #* abs(sum(log(DD[DD>0]),na.rm=TRUE)) #* abs(sum(D)) #*sd(log(val[val>0]),na.rm=TRUE)#/length(D)

#}

val_neg  <- val[neg]
val_pos  <- val[pos]

Bd_local_neg <- sum(log(val_neg[val_neg>0]),na.rm=TRUE)
Bd_local_pos <- sum(log(val_pos[val_pos>0]),na.rm=TRUE)

#Bd_local_neg <- Bd_local_neg/length(val_neg[val_neg>0])
#Bd_local_pos <- Bd_local_pos/length(val_pos[val_pos>0])


Bd_local <- abs(Bd_local_neg - Bd_local_pos) * (-1)


#Bd_local <- sum(log(val[val>0]),na.rm=TRUE) #min(Bd_local_vec, na.rm=TRUE) #sum(log(val[val>0]),na.rm=TRUE) #mean(Bd_local_vec)


# local Like from SweepFinder
#Counts1     <- table(D)
#p1          <- as.numeric(names(Counts1))
#val         <- (Counts1/sum(Counts1))^Counts1 # p^Counts  # Counts^p 
#CL[xx]      <- sum(log(val[val>0])) #prod(val[val>0])


D.names   <- as.numeric(rownames(global.mat))
#D.names   <- as.numeric(names(global.mat))

val       <- rep(NaN,length(D)-1) 
#Bd_global_vec <- rep(NaN,length(D)-1)


#for  ( yy in 1:1000){

#perm <- sample(1:length(val),length(val),replace=FALSE)

pos 	     <- rep(FALSE,length(D)-1)
neg	     <- rep(FALSE,length(D)-1)	


#perm  <- 1:(length(D)-1)

for  ( xx in perm ){

	row     <- match(D[xx],D.names)
	col     <- match(D[xx+1],D.names)

	if((D[xx]+D[xx+1])<=0){
		neg[xx] <- TRUE
	}

	if((D[xx]+D[xx+1])>=0){	
		pos[xx] <- TRUE
	}

	#xxx      <- abs(D[xx]-D[xx+1])
	#row      <- match(xxx,D.names)

	#print("D1")
	#print(D[xx])

	#print("D2")
	#print(D[xx+1])

	#print("row")
	#print(row)
	#print("col")
	#print(col)

	#valx     <- global.mat[row,col]/sum(global.mat)	

	#valx     <-  sum(global.mat[row,])/sum(global.mat)
	#valxx    <-  sum(global.mat[,col])/sum(global.mat)
	#valy     <-  global.mat[row,col]/(sum(global.mat[row,])+sum(global.mat[,col]))
		
	valz     <-    global.mat[row,col]/sum(global.mat)
       
	#valz     <-   valz * (2-abs(mean(D,na.rm=TRUE)))

	#val[xx]  <-  (valx + valxx) * valy * abs(D[xx]+D[xx+1]) #(1+abs(D[xx]+D[xx+1])) 	
	val[xx]  <-   valz #* (1-(abs(D[xx]+D[xx+1])/2)) #* (((abs(mean(D))))) #* (1-abs(D[xx])) #+D[xx+1])) #* (2-abs(D[xx]-D[xx+1])) ) 
	

	#print("value11")
	#print(value1)
	#print("value22")
	#print(value2)
	#print("valval")
	#print(val[xx])	

}

#DD        <- 0.0001 + (1-abs(D))
#Bd_global_vec[yy] <-  sum(log(val[val>0]),na.rm=TRUE) #* abs(sum(log(DD[DD>0]),na.rm=TRUE)) #* abs(sum(D)) #*sd(log(val[val>0]),na.rm=TRUE) #/length(D)
#}

val_neg  <- val[neg]
val_pos  <- val[pos]

Bd_global_neg <- sum(log(val_neg[val_neg>0]),na.rm=TRUE)
Bd_global_pos <- sum(log(val_pos[val_pos>0]),na.rm=TRUE)

#Bd_global_neg <- Bd_global_neg/length(val_neg[val_neg>0])
#Bd_global_pos <- Bd_global_pos/length(val_pos[val_pos>0])


Bd_global <- abs(Bd_global_neg - Bd_global_pos) * (-1)

#Bd_global <- sum(log(val[val>0]),na.rm=TRUE) #max(Bd_global_vec,na.rm=TRUE)

#Bd_clr     <- Bd_global
#D <- abs(D)

#P_D       <- (1-(abs(sum(D, na.rm=TRUE))/n.valid.sites))


Bd_local  <- Bd_local  #/n.valid.sites
Bd_global <- Bd_global #/n.valid.sites

Bd_clr    <- 2*( Bd_local - Bd_global )  
#Bd_clr    <- Bd_clr/sd(D)

df1       <- dim(local.mat)[1]  * dim(local.mat)[2]
df2       <- dim(global.mat)[1] * dim(global.mat)[2]

df        <- df2 - df1

P.Bd_clr  <- NaN #pchisq(Bd_clr,df=df)

SAVECLR[yy] <- Bd_clr

} # End of iteration

Bd_clr <- mean(SAVECLR, na.rm=TRUE) 

return(list(Bd_clr, P.Bd_clr))

#return(list(Bd_clr/n.valid.sites, P.Bd_clr)) #/n.valid.sites)

# CLR: global Like from SweepFinder 
#p2        <- as.numeric(names(freq.table[[xx]]))
#pids      <- match(p1,p2)
#Counts2   <- freq.table[[xx]][pids]
#val       <- (Counts2/sum(freq.table[[xx]]))^Counts1
#val       <- sum(log(val[val>0])) #prod(val[val>0])
#CLR[xx]   <- 2*( CL[xx] - val )   #2*( log(CL[xx]) - log(val) )      	
	



}
