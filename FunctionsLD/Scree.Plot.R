Scree.Plot <- function(R,main="Scree Plot",sub=NULL){
    roots <- eigen(R)$values
    x <- 1:dim(R)[1]
    plot(x,roots,type="b",col='blue',ylab="Eigenvalue",
         xlab="Component Number",main=main,sub=sub) 
    abline(h=1,lty=2,col="red")
    
}