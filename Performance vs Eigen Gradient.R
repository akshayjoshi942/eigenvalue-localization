
#Eigenvalue localization performance analysis for different Eigen Gradient
#Author: Akshay Joshi
#Description: This code implements both EL methods proposed by Cvetkovic and Chaoqian Li for random 5*5 double stochastic matrices. Can be used to test other eigengradients and variant spectrum.

#Generating random steep eigen gradients 

EValues-c(1,runif(1,0.05,0.2),runif(1,0.05,0.2),runif(1,0.05,0.2),runif(1,0.05,0.2),0)

#Extra 0 for computing elements of generated matrix. Dropped later

#Plot Parameters
myname <- c("Localization_plot_Equal_Subdominant")
mypath <- paste("C:/Users/aksha/OneDrive/Documents/NC State/Spring 2017/MA 723/Project/Plots_Eigen_Gradient/Localization_plot_Equal_subdominant_0.08", ".png", sep = "")
ColorScheme<-c("red","blue","green","pink","brown")

#Matrix generation
#Source: http://www.sciencedirect.com/science/article/pii/0024379584900181#

n=j=5

i=1
for(i in 1:5)
{
  assign(paste("alpha",i,sep=""),(EigValues[i]-EigValues[(i+1)]))
  i=i+1
}
a11<-((alpha1/n)+alpha2+alpha3+alpha4+alpha5)
a22<-(alpha1/n)+(alpha2/(n-1))+alpha3+alpha4+alpha5
a33<-(alpha1/n)+(alpha2/(n-1))+(alpha3/(n-2))+alpha4+alpha5
a44<-(alpha1/n)+(alpha2/(n-1))+(alpha3/(n-2))+(alpha4/(n-3))+alpha5
a55<-
  (alpha1/n)+(alpha2/(n-1))+(alpha3/(n-2))+(alpha4/(n-3))+(alpha5/(n-4))

a12<-a13<-a14<-a15<-a21<-a31<-a41<-a51<-alpha1/n
a23<-a24<-a25<-a32<-a42<-a52<-(alpha1/n)+(alpha2/(n-1))
a34<-a35<-a43<-a53<-(alpha1/n)+(alpha2/(n-1))+(alpha3/(n-2))
a45<-a54<-(alpha1/n)+(alpha2/(n-1))+(alpha3/(n-2))+(alpha4/(n-3))

matrix<-matrix(rbind(c(a11,a12,a13,a14,a15),
                     c(a21,a22,a23,a24,a25),
                     c(a31,a32,a33,a34,a35),
                     c(a41,a42,a43,a44,a45),
                     c(a51,a52,a53,a54,a55)), nrow=5,ncol=5)

#Dropping the extra zero
EigValues$values<-EigValues[-6]


#Summing off-diagonal elements for each row
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


#Proposed Modification by Chaoquin Li

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


#Final Table


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
legend("topleft",legend=EigValues$values,pch=24, col=ColorScheme, title="Eigenvalues")
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
legend("topleft",legend=EigValues$values,pch=24, col=ColorScheme, title="Eigenvalues")
