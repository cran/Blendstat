Blend <- function(Exp, X, Y, Conc = NULL, Effects = NULL) {
  
  # Funcao para trabalhar com problemas de otmizacao, baseada 
  # no artigo: "Kalirajan, K. P. 1990. On the estimation of a
  # regression model with fixed and random coefficients. 
  # Journal of Applied Statistics, 17: 237-244."
  # Desenvolvida por Marcelo Angelo Cirillo e 
  # Paulo Cesar Ossani em 11/2017 
  
  # Entrada:
  # Exp - Vetor com os nomes dos experimentos.
  # X - Variaveis regressoras, sem o vetor das concentracoes.
  # Y - Variavel resposta.
  # Conc - Vetor com as concentracoes dos experimentos.
  # Effects - Vetor dos efeitos das misturas em uma mistura de referencia (exemplo: centroide)
  
  # Retorna:
  # MPred  - Matriz com os valores preditos e observados.
  # MCPred - Matriz com os valores preditos por componentes.
  # MExp   - Matriz com o Design das Experiencias
  # Theta  - Vetor com as estimativas de Theta. 
 
  X   <- as.data.frame(X)
  Y   <- as.data.frame(Y)
  Exp <- as.data.frame(Exp)
  
  if (nrow(Exp) != nrow(X))
     stop("Number of lines in 'Exp' should be equal to 'X'. Verify!")
  
  if (nrow(Y) != nrow(X))
     stop("Number of lines in 'Y' should be equal to 'X'. Verify!")
  
  if (nrow(X) != length(Conc) && !is.null(Conc))
     stop("Number of lines in 'Conc' should be equal to 'X'. Verify!")
  
  if (nrow(X) != length(Effects) && !is.null(Effects))
     stop("Number of lines in 'Effects' should be equal to 'X'. Verify!")
  
  if (is.null(Effects))
     Effects <- rep(1, nrow(X))
  
  if (is.null(Conc))
     Conc <- rep(1, nrow(X))
  
  Exp.Table <- table(Exp) # tabela com as quantidade de amostras em cada experimento
  Exp.Names <- names(Exp.Table)  # nomes dos experimentos
  Num.Exp   <- length(Exp.Table) # numero de experimentos

  if ((sum(Exp.Table) / Num.Exp) != Exp.Table[[1]])
     stop("The experiments should be balanced. Verify!")
  
  Xc <- cbind(Exp, X, Conc)
  
  Y  <- cbind(Exp, Y) # acrescenta nomes dos experimentos a variavel Y
  
  ## Calculo da matris Z
  MZ <- NULL # matriz Z
  for(i in 1:Num.Exp) {
    
    MaZ <- matrix(0.0, nrow = Exp.Table[i], ncol = Num.Exp)
    MaZ[,i] <- 1
    
    MZ <- rbind(MZ, cbind(Exp.Names[i], as.data.frame(MaZ)))
    
  }
  MZ <- as.data.frame(MZ)

  ## Calculo da matris S
  MS  <- NULL # matriz S
  ncz <- ncol(MZ)
  ncx <- ncol(Xc) # numero de variaveis regressoras + 1
  for(i in 1:Num.Exp) {
    
    MaZ <- (MZ[MZ[,1] == Exp.Names[i], 2:ncz])
    
    MaX <- (Xc[Xc[,1] == Exp.Names[i], 2:ncx])
    
    MaS <- cbind(Exp.Names[i], as.data.frame(t(rbind(t(MaX),t(MaZ)))))
    
    MS  <- rbind(MS, MaS)
    
  }
  MS <- as.data.frame(MS)
  
  ## Calculo da matriz inversa e encontra as covariancias em Z
  MInv <- NULL # Matriz inversa
  MSt  <- NULL # Matriz com as covariancias de Z
  ncs  <- ncol(MS)
  # ncx  <- ncol(Xc) # numero de variaveis regressoras + 1
  for(i in 1:Num.Exp) {
    
    S <- as.matrix(MS[MS[,1] == Exp.Names[i], 2:ncs])
    
    Inv   <- ginv(t(S) %*% S)
    
    MaInv <- cbind(Exp.Names[i], as.data.frame(Inv))
    
    MaSt  <- cbind(Exp.Names[i], as.data.frame(Inv[ncx:(ncs-1), ncx:(ncs-1)])) # resultado referente aos calculo na matriz Z
    
    MInv  <- rbind(MInv, MaInv)
 
    MSt   <- rbind(MSt, MaSt)
  }
  MSt  <- as.data.frame(MSt)
  MInv <- as.data.frame(MInv)
  
  ## Calculo dos valores de theta
  MTheta <- NULL # Matriz com valores de theta
  VTheta <- NULL # Matriz com valores de theta
  nci    <- ncol(MInv)
  for(i in 1:Num.Exp) {
    
    MaI <- as.matrix(MInv[MInv[,1] == Exp.Names[i], 2:nci]) # matriz inversa de cada experimento
    MaS <- as.matrix(MS[MS[,1] == Exp.Names[i], 2:ncs]) # matriz S de cada experimento
    MaY <- as.matrix(Y[Y[,1] == Exp.Names[i], 2]) # vetor Y de cada experimento
    
    VaTheta <- MaI %*% t(MaS) %*% MaY
    
    VTheta <- cbind(VTheta, VaTheta)
    
    MaTheta <- cbind(Exp.Names[i], as.data.frame(VaTheta))
    
    MTheta <- rbind(MTheta, MaTheta)
    
  }
  MTheta <- as.data.frame(MTheta)

  
  ## Calculo dos valores de sigma
  MSigma <- NULL # Matriz com valores de sigma
  n      <- Exp.Table[[1]] # numero de elementos em cada experimento
  p      <- ncol(X) + 1 # numero de variaveis regressoras + concentracao
  q      <- 1 # numero de variaveis estocasticas
  nct    <- ncol(MTheta)
  for(i in 1:Num.Exp) {
    
    MaY  <- as.matrix(Y[Y[,1] == Exp.Names[i], 2]) # vetor Y de cada experimento
    MaS  <- as.matrix(MS[MS[,1] == Exp.Names[i], 2:ncs])  # matriz S de cada experimento
    MaTh <- as.matrix(MTheta[MTheta[,1] == Exp.Names[i], 2:nct]) # matriz theta de cada experimento
    
    Ma      <- (MaY - MaS%*%MaTh) 
    
    MaSigma <- cbind(Exp.Names[i], as.data.frame((t(Ma) %*% Ma) / (n - p - q)))
    
    MSigma  <- rbind(MSigma, MaSigma)
    
  }
  MSigma <- as.data.frame(MSigma)

  ## Calculo dos valores de delta
  T1 <- cov(VTheta)
  T2 <- 0
  ncst <- ncol(MSt)
  for(i in 1:Num.Exp) {

    MaSt <- as.matrix(MSt[MSt[,1] == Exp.Names[i], 2:ncst]) # matriz theta de cada experimento
    
    T2 <- T2 + MSigma[i,2] * MaSt
  }
  
  T2 <- T2 / Num.Exp 
  
  dif <- T1 - T2 
  
  dvs <- svd(dif)
  
  if (length(dvs$d) > 1) {
     MDelta <- diag(dvs$d)
  } else {
     MDelta <- dvs$d
  }
  
  ## Calculo da matris V
  MV   <- NULL # matriz V
  t    <- NULL
  auxV <- matrix (0, nrow = n, ncol = n)
  for(i in 1:Num.Exp) {
    
    MaZ <- as.matrix(MZ[MZ[,1] == Exp.Names[i], 2:ncz])
 
    V <- MaZ %*% MDelta %*% t(MaZ) + MSigma[i,2] * diag(1, n)
    
    for(j in 1:Num.Exp) {
      
      if (j == i) { aux = V }
      else aux = auxV
      
      t <- cbind(t, aux)
    }

    MV <- rbind(MV, t)
    
    t <- NULL
  }
  MVInv <- ginv(MV)

  ## Estimativa de Theta
  NS <- NULL
  ncs <- ncol(MS)
  for(i in 1:Num.Exp) {
    
    MaS <- as.matrix(MS[MS[,1] == Exp.Names[i], 2:ncs]) # matriz S de cada experimento

    NS <- rbind(NS, MaS)
  }

  MInvTheta <- ginv(t(NS) %*% MVInv %*% NS)

  EstTheta <- MInvTheta %*% t(NS) %*% MVInv %*% Y[,2]
  colnames(EstTheta) <- "Theta"
  rownames(EstTheta) <- c(colnames(Xc[,2:ncol(Xc)]), paste("Exp", Exp.Names))

  ## Calculo dos Valores preditos
  MU <- as.matrix(cbind(Xc[,2:ncol(Xc)], MZ[,2:ncol(MZ)]))
  VlrPred <- MU %*% EstTheta
  
  ## Matriz com os valores Preditos
  MPred <- cbind(Exp, Effects, Conc, VlrPred, Y[,2])
  colnames(MPred) <- c("Experiments", "Effects", "Concentrations", "Predicted values", "Observed values")
  MPred <- as.data.frame(MPred)

  ## Matriz com os valores Preditos por Componentes
  MCPred <- NULL
  X <- as.matrix(X)
  nreg <- ncol(X) # numero de variaveis regressoras
  for(i in 1:nreg) { 
    MCPred <- cbind(MCPred, X[,i] * EstTheta[i] + EstTheta[nreg + 1] * Conc) # variaveis regressoras + concentracoes
    for(j in 1:Num.Exp) {
      MCPred[,i] <- MCPred[,i] + EstTheta[nreg + 1 + j] * MZ[,j+1] # + efeitos da matriz Z
    }
  }
  MCPred <- cbind(Exp, Effects, MCPred)  
  colnames(MCPred) <- c("Experiments", "Effects", colnames(X))
  MCPred <- as.data.frame(MCPred)

  ## Matriz de experimentos
  MExp <- as.matrix(MZ[,2:ncol(MZ)])
  colnames(MExp) <- paste("Exp.", 1:ncol(MExp))
  
  Lista <- list(MPred = MPred, MCPred = MCPred, Theta = EstTheta, MExp = MExp)

  return(Lista)
}
