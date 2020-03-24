Blend <- function(exp, X, Y, conc = NULL, effects = NULL) {
  
  # Funcao para trabalhar com problemas de otmizacao, baseada 
  # no artigo: "Kalirajan, K. P. 1990. On the estimation of a
  # regression model with fixed and random coefficients. 
  # Journal of Applied Statistics, 17: 237-244."
  # Desenvolvida por Marcelo Angelo Cirillo e 
  # Paulo Cesar Ossani em 11/2017 
  
  # Entrada:
  # exp - Vetor com os nomes dos experimentos.
  # X - Variaveis regressoras, sem o vetor das concentracoes.
  # Y - Variavel resposta.
  # conc - Vetor com as concentracoes dos experimentos.
  # effects - Vetor dos efeitos das misturas em uma mistura de referencia (exemplo: centroide)
  
  # Retorna:
  # MPred  - Matriz com os valores preditos e observados.
  # MCPred - Matriz com os valores preditos por componentes.
  # Mexp   - Matriz com o Design das experiencias
  # theta  - Vetor com as estimativas de theta. 
 
  X   <- as.data.frame(X)
  Y   <- as.data.frame(Y)
  exp <- as.data.frame(exp)
  
  if (nrow(exp) != nrow(X))
     stop("Number of lines in 'exp' should be equal to 'X'. Verify!")
  
  if (nrow(Y) != nrow(X))
     stop("Number of lines in 'Y' should be equal to 'X'. Verify!")
  
  if (nrow(X) != length(conc) && !is.null(conc))
     stop("Number of lines in 'conc' should be equal to 'X'. Verify!")
  
  if (nrow(X) != length(effects) && !is.null(effects))
     stop("Number of lines in 'effects' should be equal to 'X'. Verify!")
  
  if (is.null(effects))
     effects <- rep(1, nrow(X))
  
  if (is.null(conc))
     conc <- rep(1, nrow(X))
  
  exp.Table <- table(exp) # tabela com as quantidade de amostras em cada experimento
  exp.Names <- names(exp.Table)  # nomes dos experimentos
  num.exp   <- length(exp.Table) # numero de experimentos

  if ((sum(exp.Table) / num.exp) != exp.Table[[1]])
     stop("The experiments should be balanced. Verify!")
  
  Xc <- cbind(exp, X, conc)
  
  Y  <- cbind(exp, Y) # acrescenta nomes dos experimentos a variavel Y
  
  ## Calculo da matris Z
  MZ <- NULL # matriz Z
  for(i in 1:num.exp) {
    
    MaZ <- matrix(0.0, nrow = exp.Table[i], ncol = num.exp)
    MaZ[,i] <- 1
    
    MZ <- rbind(MZ, cbind(exp.Names[i], as.data.frame(MaZ)))
    
  }
  MZ <- as.data.frame(MZ)

  ## Calculo da matris S
  MS  <- NULL # matriz S
  ncz <- ncol(MZ)
  ncx <- ncol(Xc) # numero de variaveis regressoras + 1
  for(i in 1:num.exp) {
    
    MaZ <- (MZ[MZ[,1] == exp.Names[i], 2:ncz])
    
    MaX <- (Xc[Xc[,1] == exp.Names[i], 2:ncx])
    
    MaS <- cbind(exp.Names[i], as.data.frame(t(rbind(t(MaX),t(MaZ)))))
    
    MS  <- rbind(MS, MaS)
    
  }
  MS <- as.data.frame(MS)
  
  ## Calculo da matriz inversa e encontra as covariancias em Z
  MInv <- NULL # Matriz inversa
  MSt  <- NULL # Matriz com as covariancias de Z
  ncs  <- ncol(MS)
  # ncx  <- ncol(Xc) # numero de variaveis regressoras + 1
  for(i in 1:num.exp) {
    
    S <- as.matrix(MS[MS[,1] == exp.Names[i], 2:ncs])
    
    Inv   <- ginv(t(S) %*% S)
    
    MaInv <- cbind(exp.Names[i], as.data.frame(Inv))
    
    MaSt  <- cbind(exp.Names[i], as.data.frame(Inv[ncx:(ncs-1), ncx:(ncs-1)])) # resultado referente aos calculo na matriz Z
    
    MInv  <- rbind(MInv, MaInv)
 
    MSt   <- rbind(MSt, MaSt)
  }
  MSt  <- as.data.frame(MSt)
  MInv <- as.data.frame(MInv)
  
  ## Calculo dos valores de theta
  Mtheta <- NULL # Matriz com valores de theta
  Vtheta <- NULL # Matriz com valores de theta
  nci    <- ncol(MInv)
  for(i in 1:num.exp) {
    
    MaI <- as.matrix(MInv[MInv[,1] == exp.Names[i], 2:nci]) # matriz inversa de cada experimento
    MaS <- as.matrix(MS[MS[,1] == exp.Names[i], 2:ncs]) # matriz S de cada experimento
    MaY <- as.matrix(Y[Y[,1] == exp.Names[i], 2]) # vetor Y de cada experimento
    
    Vatheta <- MaI %*% t(MaS) %*% MaY
    
    Vtheta <- cbind(Vtheta, Vatheta)
    
    Matheta <- cbind(exp.Names[i], as.data.frame(Vatheta))
    
    Mtheta <- rbind(Mtheta, Matheta)
    
  }
  Mtheta <- as.data.frame(Mtheta)

  
  ## Calculo dos valores de sigma
  MSigma <- NULL # Matriz com valores de sigma
  n      <- exp.Table[[1]] # numero de elementos em cada experimento
  p      <- ncol(X) + 1 # numero de variaveis regressoras + concentracao
  q      <- 1 # numero de variaveis estocasticas
  nct    <- ncol(Mtheta)
  for(i in 1:num.exp) {
    
    MaY  <- as.matrix(Y[Y[,1] == exp.Names[i], 2]) # vetor Y de cada experimento
    MaS  <- as.matrix(MS[MS[,1] == exp.Names[i], 2:ncs])  # matriz S de cada experimento
    MaTh <- as.matrix(Mtheta[Mtheta[,1] == exp.Names[i], 2:nct]) # matriz theta de cada experimento
    
    Ma      <- (MaY - MaS%*%MaTh) 
    
    MaSigma <- cbind(exp.Names[i], as.data.frame((t(Ma) %*% Ma) / (n - p - q)))
    
    MSigma  <- rbind(MSigma, MaSigma)
    
  }
  MSigma <- as.data.frame(MSigma)

  ## Calculo dos valores de delta
  T1 <- cov(Vtheta)
  T2 <- 0
  ncst <- ncol(MSt)
  for(i in 1:num.exp) {

    MaSt <- as.matrix(MSt[MSt[,1] == exp.Names[i], 2:ncst]) # matriz theta de cada experimento
    
    T2 <- T2 + MSigma[i,2] * MaSt
  }
  
  T2 <- T2 / num.exp 
  
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
  for(i in 1:num.exp) {
    
    MaZ <- as.matrix(MZ[MZ[,1] == exp.Names[i], 2:ncz])
 
    V <- MaZ %*% MDelta %*% t(MaZ) + MSigma[i,2] * diag(1, n)
    
    for(j in 1:num.exp) {
      
      if (j == i) { aux = V }
      else aux <- auxV
      
      t <- cbind(t, aux)
    }

    MV <- rbind(MV, t)
    
    t <- NULL
  }
  MVInv <- ginv(MV)

  ## Estimativa de theta
  NS  <- NULL
  ncs <- ncol(MS)
  for(i in 1:num.exp) {
    
    MaS <- as.matrix(MS[MS[,1] == exp.Names[i], 2:ncs]) # matriz S de cada experimento

    NS <- rbind(NS, MaS)
  }

  MInvtheta <- ginv(t(NS) %*% MVInv %*% NS)

  Esttheta <- MInvtheta %*% t(NS) %*% MVInv %*% Y[,2]
  colnames(Esttheta) <- "theta"
  rownames(Esttheta) <- c(colnames(Xc[,2:ncol(Xc)]), paste("exp", exp.Names))

  ## Calculo dos Valores preditos
  MU <- as.matrix(cbind(Xc[,2:ncol(Xc)], MZ[,2:ncol(MZ)]))
  VlrPred <- MU %*% Esttheta
  
  ## Matriz com os valores Preditos
  MPred <- cbind(exp, effects, conc, VlrPred, Y[,2])
  colnames(MPred) <- c("Experiments", "Effects", "Concentrations", "Predicted values", "Observed values")
  MPred <- as.data.frame(MPred)

  ## Matriz com os valores Preditos por Componentes
  MCPred <- NULL
  X <- as.matrix(X)
  nreg <- ncol(X) # numero de variaveis regressoras
  for(i in 1:nreg) { 
    MCPred <- cbind(MCPred, X[,i] * Esttheta[i] + Esttheta[nreg + 1] * conc) # variaveis regressoras + concentracoes
    for(j in 1:num.exp) {
      MCPred[,i] <- MCPred[,i] + Esttheta[nreg + 1 + j] * MZ[,j+1] # + efeitos da matriz Z
    }
  }
  MCPred <- cbind(exp, effects, MCPred)  
  colnames(MCPred) <- c("experiments", "effects", colnames(X))
  MCPred <- as.data.frame(MCPred)

  ## Matriz de experimentos
  Mexp <- as.matrix(MZ[,2:ncol(MZ)])
  colnames(Mexp) <- paste("exp.", 1:ncol(Mexp))
  
  Lista <- list(MPred = MPred, MCPred = MCPred, theta = Esttheta, Mexp = Mexp)

  return(Lista)
}
