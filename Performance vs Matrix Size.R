#Eigenvalue localization analysis for increasing matrix sizes
#Author: Akshay Joshi

#This code creates random row stochastic matrices from size 5 to 100 in increments of 5. 
#Runs both the Cvetkovic and Chaoqian's method for eigenvalue localization and saves the final plot in specified folder. 

#Row Creation
i=1
j=1
for(j in 5*(1:20))
{
  for (i in 1:j)
  {
    x<-runif(j, min=0)
    assign(paste("row", i, sep = ""), x/sum(x))
    i=i+1
  }
  
  #Build Matrix
  i=1
  MatrixBuild=rbind(row1,row2)
  for (i in 3:j)
  {
    
    assign(paste("MatrixBuild"),rbind(MatrixBuild,get(paste("row", i, sep = ""))))
    i=i+1
  }
  
  #Matrix Set up
  matrix<-matrix(MatrixBuild, nrow = j, ncol = j)
  matrix
  EigValues<-as.data.frame(eigen(matrix,only.values = FALSE))
  EigValues$values[]
  
  #Final Plot Parameters
  ColorScheme<-sample(colors()[-1],j)
  mypath <- paste("C:/Users/aksha/OneDrive/Documents/NC State/Spring 2017/MA 723/Project/Plots_Matrix_Size/Final_plot_Stochastic_Matrix_Size_", j, ".png", sep = "")
  myname <- paste("Localization_plot_Stochastic_Matrix_Size_", j, sep = "")
  
  i=1
  for (i in 1:j)
  {
    assign(paste("r", i, sep = ""), 1-matrix[i,i])
    i=i+1
  }
  
  #Cvetkovic's Method
  #S i values
  i=1
  for (i in 1:j)
  {
    
    { 
      assign(paste("s", i, sep = ""),min(matrix[-i,i]))
      i=i+1
    }
  }
  
  #Gamma Build
  
  i=1
  for (i in 1:j)
  {
    
    { 
      assign(paste("Gamma", i, sep = ""),matrix[i,i]-get(paste("s", i, sep = ""  )))
      i=i+1
    }
  }
  
  #Gamma Value
  i=1
  GammaValue<-max(Gamma1,Gamma2)
  for (i in 1:(j))
  {
    
    { 
      assign("GammaValue",max(get(paste("GammaValue")),get(paste("Gamma", i, sep = ""))))
      i=i+1
      print(GammaValue)
    }
  }
  
  
  #r(matrix)
  rmatrix<-1-sum(diag(matrix))+((j-1)*GammaValue)
  
  #Improved Method
  #S i values
  i=1
  for (i in 1:j)
  {
    
    { 
      assign(paste("Ss", i, sep = ""),max(matrix[i,-i]))
      i=i+1
    }
  }
  
  
  i=1
  for (i in 1:j)
  {
    
    { 
      assign(paste("GammaTilde", i, sep = ""),matrix[i,i]-get(paste("Ss", i, sep = ""  )))
      i=i+1
    }
  }
  
  
  i=1
  GammaTildeValue<-max(GammaTilde1,GammaTilde2)
  for (i in 1:(j))
  {
    
    { 
      assign("GammaTildeValue",max(get(paste("GammaTildeValue")),get(paste("GammaTilde", i, sep = ""))))
      i=i+1
      print(GammaTildeValue)
    }
  }
  GammaTildeValue
  
  #r(matrix)
  rTildematrix<-sum(diag(matrix))+((j-1)*GammaTildeValue)-1
  
  #Checking if the eigenvalues lie inside the Improved Method Circle
  
  Improved_Method<-rep(0,times=j)
  Eigenvalue<-rep(0,times=j)
  Results_Table<-as.data.frame(cbind(Eigenvalue,Improved_Method))
  Results_Table$Improved_Distance<-rep(0,times=j)
  
  
  for (i in 1:j)
  {
    Results_Table$Eigenvalue[i]<-EigValues$value[i]
    Results_Table$Improved_Distance[i]<-(sqrt(((Re(EigValues$value[i])-Re(GammaTildeValue))^2)+ ((Im(EigValues$value[i])-Im(GammaTildeValue))^2)))
    Results_Table$Improved_Method[i]<-Results_Table$Improved_Distance[i]<rTildematrix
  }
  
  Results_Table$Improved_Method<-as.numeric(Results_Table$Improved_Method)
  Results_Table
  
  #Checking if the eigenvalues lie inside the Cvetkovic's Circle
  
  Results_Table$Cvetkovic_Method<-rep(0,times=j)
  Results_Table$Cvetkovic_Distance<-rep(0,times=j)
  
  for (i in 1:j)
  {
    Results_Table$Cvetkovic_Distance[i]<-(sqrt(((Re(EigValues$value[i])-Re(GammaValue))^2)+ ((Im(EigValues$value[i])-Im(GammaValue))^2)))
    Results_Table$Cvetkovic_Method[i]<-Results_Table$Cvetkovic_Distance[i]<rmatrix
  }
  
  Results_Table
  Results_Table<-Results_Table[which(Re(Results_Table$Eigenvalue) < 1),]
  
  #Final Plot
  library("plotrix")
  
  png(file=mypath)
  
  plot(0,type='n',main=myname,xlim=c(-2,2), ylim = c(-2,2), xlab = "Real", ylab = "Imaginary")
  
  i=1
  for (i in 1:j)
  {
    draw.circle(matrix[i,i],0,get(paste("r", i, sep = ""  )),border =ColorScheme[i],lwd=2)
    points(Re(EigValues$values[i]),Im(EigValues$values[i]),col=ColorScheme[i],pch=24)
    i=i+1
  }
  draw.circle(GammaValue,0,rmatrix,lty=3,lwd=3)
  draw.circle(GammaTildeValue,0,rTildematrix,lty=5,lwd=3, border = "grey")
  dev.off()
  
  plot(0,type='n',main=myname,xlim=c(-2,2), ylim = c(-2,2), xlab = "Real", ylab = "Imaginary")
  
  i=1
  for (i in 1:j)
  {
    draw.circle(matrix[i,i],0,get(paste("r", i, sep = ""  )),border =ColorScheme[i],lwd=2)
    points(Re(EigValues$values[i]),Im(EigValues$values[i]),col=ColorScheme[i],pch=24)
    i=i+1
  }
  draw.circle(GammaValue,0,rmatrix,lty=3,lwd=3)
  draw.circle(GammaTildeValue,0,rTildematrix,lty=5,lwd=3, border = "grey")
  
  j=j+1
  rm(list=ls())
}


