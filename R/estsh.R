###############################################################
######################code for R package#######################
###############################################################

#' Generate an Latin square (LS)
#'
#' @param d an integer, the run size of the resulting LS.
#'
#' @return an LS with d rows and d columns.
#' @export
#'
#' @examples
#' onels(5)
onels<-function(d){
  sq <- matrix(nrow=d, ncol=d)
  sq[1,]<-c(1,sample(d-1)+1);
  for (r in 2:d) {
    sq[r,]<-c(sq[r-1,-1],sq[r-1,1])
  }
  pcol<-sample(d);prow<-c(1,sample(c(2:d)));
  sq<-sq[,pcol]; sq<-sq[prow,];
  return(sq)
}


#' Estimating the Shapley value based on Latin square (LS)
#'
#' @param d an integer, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param ... other parameters used in val(sets,...).
#'
#' @return a vector including estimated Shapley values of all players based on LS.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=8,ncol=8)
#' temp_adjacent[1,6:8]<-1;temp_adjacent[2,7]<-1;temp_adjacent[c(4,6,7),8]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.shls(8,56,temp_val,temp_adjacent)

est.shls<-function(d,n,val,...){
  sh<-rep(0,d); k<-n/d;
  if(k%%1!=0){stop('Error: n should be a multiple of d')}
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],...)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


#' Determine whether an integer is a prime
#'
#' @param x the integer to be determined.
#'
#' @return the result: TRUE (x is a prime) or FALSE (x is not a prime).
#' @export
#'
#' @examples
#' is.prime(7)
#' is.prime(8)
is.prime<-function(x){
  I<-1
  if(!x%in%c(2,3)){
    for(i in 2:floor(x/2)){
      if(x%%i==0){I<-0;break;}
    }
  }
  if(I==1){return(TRUE)}
  else{return(FALSE)}
}


#' Generate a component orthogonal array (COA) with a prime d
#'
#' @param d a prime, the column of the resulting COA.
#'
#' @return a COA with d(d-1) rows and d columns.
#' @export
#'
#' @examples
#' onecoa.prime(5)
onecoa.prime<-function(d){
  if(!is.prime(d)) {stop('Error: d should be a prime')}
  firstline<-c(0,sample(d-1)); cz<-matrix(0,d-1,d); D<-matrix(0,d*(d-1),d)
  for(i in 1:(d-1)){
    cz[i,]<-(firstline*i)%%d
  }
  for(j in 0:(d-1)){
    D[(j*(d-1)+1):((j+1)*(d-1)),]<-(cz+j)%%d
  }
  D<-D+1
  return(D)
}

#' Estimating the Shapley value based on component orthogonal array (COA) with a prime d
#'
#' @param d a prime, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param ... other parameters used in val(sets,...).
#'
#' @return a vector including estimated Shapley values of all players based on COA.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=5,ncol=5)
#' temp_adjacent[1,c(2,3,5)]<-1;temp_adjacent[2,4]<-1;temp_adjacent[3,5]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.shcoa.prime(5,20,temp_val,temp_adjacent)
est.shcoa.prime<-function(d,n,val,...){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  if(k%%1!=0){stop('Error: n should be a multiple of d(d-1)')}
  for(t in 1:k){
    coat<-onecoa.prime(d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],...)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


#' Polynomial division
#'
#' @param f1 a vector represents the coefficients of the dividend polynomial. For example, if the dividend is x^5+2x^3+1, then f1=c(1,0,2,0,0,1).
#' @param f2 a vector represents the coefficients of the dividend polynomial. For example, if the divisor is x^4+2, then f2=c(1,0,0,0,2).
#'
#' @return a vector represents the coefficients of the resulting polynomial. For example, the result c(2,0,-2,1) represents 2x^3-2x+1.
#' @export
#'
#' @examples
#' poly.div(c(1,0,2,0,0,1),c(1,0,0,0,2))
poly.div<-function(f1,f2){
  d1<-length(f1)-1; d2<-length(f2)-1;
  if(d1<d2){f1<-c(rep(0,d2-d1),f1); d1<-d2}
  dd<-d1-d2;
  divisend<-f1; divisor<-f2;
  b<-rep(0,dd+1)
  for(i in 1:(dd+1)){
    b[i]<-divisend[1]/divisor[1]
    for(j in 1:(d2+1)){
      divisend[j]<-divisend[j]-b[i]*divisor[j]
    }
    divisend<-divisend[-1]
  }
  r<-divisend
  return(r)
}

#' Polynomial division defined on GF(s) with a prime s
#'
#' @param f1 a vector represents the coefficients of the dividend polynomial. For example, if the dividend is x^5+2x^3+1, then f1=c(1,0,2,0,0,1).
#' @param f2 a vector represents the coefficients of the dividend polynomial. For example, if the divisor is x^4+2, then f2=c(1,0,0,0,2).
#' @param s a prime, the order of the Galois filed (GF).
#'
#' @return a vector represents the coefficients of the resulting polynomial. For example, the result c(2,0,1,1) represents 2x^3+x+1.
#' @export
#'
#' @examples
#' gfpoly.div(c(1,0,2,0,0,1),c(1,0,0,0,2),3)
gfpoly.div<-function(f1,f2,s){
  r<-poly.div(f1,f2)
  gfr<-r%%s
  return(gfr)
}


#' Polynomial mutiplication defined on GF(s) with a prime s
#'
#' @param f1 a vector represents the coefficients of the first multiplier polynomial. For example, if the dividend is x^5+2x^3+1, then f1=c(1,0,2,0,0,1).
#' @param f2 a vector represents the coefficients of the second multiplier polynomial. For example, if the divisor is x^4+2, then f2=c(1,0,0,0,2).
#' @param s a prime, the order of the Galois filed (GF).
#'
#' @return a vector represents the coefficients of the resulting polynomial. For example, the result c(1,0,2,0,2,1,1,0,0,2) represents x^9+2x^7+2x^5+x^4+x^3+2.
#' @export
#'
#' @examples
#' gfpoly.multi(c(1,0,2,0,0,1),c(1,0,0,0,2),3)
gfpoly.multi<-function(f1,f2,s){ #f1,f2 are multipliers
  d1<-length(f1);d2<-length(f2);hd<-d1+d2-1;
  t<-matrix(0,d1,hd)
  for(i in 1:d1){ # each coefficient of f1 multiply d_2
    t[i,i:(i+d2-1)]<-f1[i]*f2 #the second coefficient is one place behind the first
  }
  multi<-apply(t,2,sum)
  gft<-multi%%s #GF(s)
  return(gft)
}


#' Polynomial additive defined on GF(s) with a prime s
#'
#' @param f1 a vector represents the coefficients of the first addend polynomial. For example, if the dividend is x^5+2x^3+1, then f1=c(1,0,2,0,0,1).
#' @param f2 a vector represents the coefficients of the second addend polynomial. For example, if the divisor is x^4+2, then f2=c(1,0,0,0,2).
#' @param s a prime, the order of the Galois filed (GF).
#'
#' @return a vector represents the coefficients of the resulting polynomial. For example, the result c(1, 1, 2, 0, 0, 0) represents x^5+x^4+2x^3.
#' @export
#'
#' @examples
#' gfpoly.add(c(1,0,2,0,0,1),c(1,0,0,0,2),3)
gfpoly.add<-function(f1,f2,s){
  l1<-length(f1); l2<-length(f2);
  if(l1<l2){f1<-c(rep(0,l2-l1),f1);}
  else if(l1>l2){f2<-c(rep(0,l1-l2),f2);}
  f<-(f1+f2)%%s
  return(f)
}

#' Generate a component orthogonal array (COA) with a prime power d
#'
#' @param d a power of prime p, the column of the resulting COA.
#' @param p a prime, the bottom number of d.
#' @param f_d a vector represents the coefficients of primative polynomial on GF(d). For example the primative polynomial on GF(3^2) is x^2+x+2, then let f_d=c(1,1,2).
#'
#' @return a COA with d(d-1) rows and d columns.
#' @export
#'
#' @examples
#' onecoa(9,3,c(1,1,2))
onecoa<-function(d,p,f_d){
  if(!is.prime(p)){stop('Error: p should be a prime.')}
  firstr<-c(0,sample(d-1));
  r<-round(log(d,p))
  #multiplication table and addition table on GF(d)
  if(r==1){
    M_d<-matrix(0,p,p); A_d<-matrix(0,p,p);
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        M_d[i,j]<-((i-1)*(j-1))%%p; M_d[j,i]<-M_d[i,j];
        A_d[i,j]<-(i+j-2)%%p; A_d[j,i]<-A_d[i,j];
      }
      M_d[i,i]<-(i-1)^2%%p; A_d[i,i]<-(2*(i-1))%%p;
    }
    M_d[p,p]<-(p-1)^2%%p; A_d[p,p]<-(2*(p-1))%%p;
  }
  else{
    entry<-gtools::permutations(p,r,0:(p-1),repeats=TRUE); cases<-gtools::permutations(d,2,repeats= TRUE);
    M1<-matrix(nrow=d,ncol=d*r)
    for(i in 1:(d^2)){
      M1[cases[i,1],((cases[i,2]-1)*r+1):(cases[i,2]*r)]<-gfpoly.div(gfpoly.multi(entry[cases[i,1],],entry[cases[i,2],],p),f_d,p)
    }
    A1<-matrix(nrow=d,ncol=d*r)
    for(i in 1:(d^2)){
      A1[cases[i,1],((cases[i,2]-1)*r+1):(cases[i,2]*r)]<-gfpoly.add(entry[cases[i,1],],entry[cases[i,2],],p)
    }
    M_d<-matrix(nrow=d,ncol=d); A_d<-matrix(nrow=d,ncol=d)
    for(i in 1:d){
      for(j in 1:d){
        I<-0
        for(l in 1:d){
          if(round(sum(abs(M1[i,((j-1)*r+1):(j*r)]-entry[l,])))==0) {M_d[i,j]<-l-1; I<-I+1;}
          if(round(sum(abs(A1[i,((j-1)*r+1):(j*r)]-entry[l,])))==0) {A_d[i,j]<-l-1; I<-I+1;}
          if(I==2) break;
        }
      }
    }
  }
  #construction
  nr<-d*(d-1); coa<-matrix(nrow=nr,ncol=d);
  for(i in 1:d){
    for(j in 1:(d-1)){
      for(k in 1:d)
        coa[(i-1)*(d-1)+j,k]<-A_d[i,M_d[j+1,firstr[k]+1]+1]
    }
  }
  return(coa+1)
}


#' Estimating the Shapley value based on component orthogonal array (COA) with a prime power d
#'
#' @param d a power of prime p, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param p a prime, the bottom number of d.
#' @param f_d a vector represents the coefficients of primative polynomial on GF(d). For example the primative polynomial on GF(3^2) is x^2+x+2, then let f_d=c(1,1,2).
#' @param ... other parameters used in val(sets,...).
#'
#' @return a vector including estimated Shapley values of all players based on COA.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=8,ncol=8)
#' temp_adjacent[1,6:8]<-1;temp_adjacent[2,7]<-1;temp_adjacent[c(4,6,7),8]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.shcoa(8,112,temp_val,2,c(1,0,1,1),temp_adjacent)
est.shcoa<-function(d,n,val,p,f_d,...){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  if(k%%1!=0){stop('Error: n should be a multiple of d(d-1)')}
  if((length(f_d)!=log(d,p)+1)|(max(f_d)>=p)){stop('Error: f_d is not a correct primitive polynomial')}
  for(t in 1:k){
    coat<-onecoa(d,p,f_d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],...)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

#' Estimating the Shapley value based on simple random sampling (SRS)
#'
#' @param d an integer, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param ... other parameters used in val(sets,...).
#'
#' @return a vector including estimated Shapley values of all players based on SRS.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=8,ncol=8)
#' temp_adjacent[1,6:8]<-1;temp_adjacent[2,7]<-1;temp_adjacent[c(4,6,7),8]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.shsrs(8,112,temp_val,temp_adjacent)
est.shsrs<-function(d,n,val,...){
  sh<-rep(0,d)
  for(l in 1:n){
    perml<-sample(d); preC<-0;
    for(i in 1:d){
      delta<-val(perml[1:i],...)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

#' Generate the structured samples of simple random samples
#'
#' @param permatrix a matrix, each row is a permutation.
#' @param jcom an integer, represents the target component. Hope that the component jcom appears the same number of at each position.
#' @param d the number of components.
#'
#' @return a matrix represents the structured samples.
#' @export
#'
#' @examples
#' temp_samples<-matrix(nrow=10,ncol=5)
#' for(i in 1:10){temp_samples[i,]<-sample(1:5,5)}
#' structed.perm(temp_samples,3,5)
structed.perm<-function(permatrix,jcom,d){
  nr<-nrow(permatrix); t<-nr/d
  jpermatrix<-matrix(0,nrow=nr,ncol=d)
  for(i in 1:nr){
    tt<-ceiling(i/t); perm<-permatrix[i,]
    id<-which(perm==jcom)
    tcom<-perm[tt]
    perm[tt]<-jcom
    perm[id]<-tcom
    jpermatrix[i,]<-perm
  }
  return(jpermatrix)
}

#' Estimating the Shapley value based on structured simple random sampling (StrRS)
#'
#' @param d an integer, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param ... other parameters used in val(sets,...).
#'
#' @return a vector including estimated Shapley values of all players based on StrRS.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=8,ncol=8)
#' temp_adjacent[1,6:8]<-1;temp_adjacent[2,7]<-1;temp_adjacent[c(4,6,7),8]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.shstrrs(8,112,temp_val,temp_adjacent)
est.shstrrs<-function(d,n,val,...){
  sh<-rep(0,d); groupsize<-n/d
  srsperm<-matrix(0,nrow=n,ncol=d)
  for(i in 1:n){
    srsperm[i,]<-sample(d)
  }
  for(j in 1:d){
    jperm<-structed.perm(srsperm,j,d)
    for(i in (groupsize+1):n){
      loc<-ceiling(i/groupsize)
      sh[j]<-sh[j]+val(jperm[i,1:loc],...)-val(jperm[i,1:(loc-1)],...)
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


#' The main algorithm for estimating the Shapley value
#'
#' @param method the method used for estimating,'SRS' meas simple random sampling, 'StrRS' means structured simple random sampling, 'LS' means Latin square and 'COA' means component orthogonal array.
#' @param d an integer, the number of players.
#' @param n an integer, the sample size.
#' @param val the predefined value function.
#' @param ... other parameters used in val(sets,...).
#' @param p a prime, the bottom number of d.
#' @param f_d a vector represents the coefficients of primative polynomial on GF(d). For example the primative polynomial on GF(3^2) is x^2+x+2, then let f_d=c(1,1,2).
#'
#' @return a vector including estimated Shapley values of all players.
#' @export
#'
#' @examples
#' temp_adjacent<-matrix(0,nrow=8,ncol=8)
#' temp_adjacent[1,6:8]<-1;temp_adjacent[2,7]<-1;temp_adjacent[c(4,6,7),8]<-1;
#' temp_adjacent<-temp_adjacent+t(temp_adjacent)
#' temp_val<-function(sets,adjacent){
#'   if(length(sets)==1) val<-0
#'   else{
#'     subadjacent<-adjacent[sets,sets]
#'     nsets<-length(sets)
#'     A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
#'     for(l in 1:(nsets-1)){
#'       A<-A%*%subadjacent
#'       B<-B+A
#'     }
#'     val<-ifelse(sum(B==0)>nsets,0,1)
#'   }
#'   return(val)
#' }
#' est.sh('SRS',8,112,temp_val,temp_adjacent)
#' est.sh('StrRS',8,112,temp_val,temp_adjacent)
#' est.sh('LS',8,112,temp_val,temp_adjacent)
#' est.sh('COA',8,112,temp_val,temp_adjacent,p=2,f_d=c(1,0,1,1))
est.sh<-function(method,d,n,val,...,p=NA,f_d=NA){
  if(method=="LS"){sh<-est.shls(d,n,val,...)}
  else if(method=="COA"){
    if(is.prime(d)){sh<-est.shcoa.prime(d,n,val,...)}
    else if(!is.na(p)){
      if(log(d,p)%%1==0){sh<-est.shcoa(d,n,val,p,f_d,...)}
      else{stop('Error: d is not a power of p')}
    }
    else{stop('Error: d is not a prime power, or p and f_d are not given')}
  }
  else if(method=="SRS"){sh<-est.shsrs(d,n,val,...)}
  else if(method=="StrRS"){sh<-est.shstrrs(d,n,val,...)}
  else{
    stop('Error: method should be LS, COA, SRS or StrRS')
  }
  return(sh)
}

