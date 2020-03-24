Plot.Blend <- function(BL, titles = c(NA,NA), posleg = 2, xlabel = NA, 
                       ylabel = NA, boxleg = FALSE, color = TRUE, 
                       expcolor = NA, casc = TRUE) {
  
  # Rotina para plotar graficos de Blendstate,
  # Desenvolvida por Marcelo Angelo Cirillo e
  # Paulo Cesar Ossani em 11/2017 
  
  # Entrada:
  # BL     - Dados da funcao Blend.
  # titles - Titulos para o grafico dos efeitos das concentracoes e componentes. Se nao for definido assume texto padrao.
  # posleg - 1 para legenda a esquerda,
  #          2 para legenda a direita (default),
  #          3 para legenda acima,
  #          4 para legenda abaixo.
  # xlabel - Nomeia o eixo X, se nao definido retorna padrao.
  # ylabel - Nomeia o eixo Y, se nao definido retorna padrao.
  # boxleg - Coloca moldura na legenda (default = TRUE).
  # color  - Graficos coloridos (default = TRUE).
  # expcolor - Vetor com as cores dos experimentos.
  # casc   - Efeito cascata na apresentacao dos graficos (default = TRUE).

  # Retorna:
  # Varios graficos.
  
  if (!is.numeric(posleg) || posleg < 1 || posleg > 4 || (floor(posleg)-posleg) != 0)
     stop("Input to set the position of the legend 'posleg' is incorrect, should be an integer number between [1,4]. Verify!")
  
  if (!is.logical(boxleg)) 
     stop("'boxleg' input is incorrect, it should be TRUE or FALSE. Verify!")
  
  if (!is.character(xlabel) && !is.na(xlabel))
     stop("'xlabel' input is incorrect, it should be of type character or string. Verify!")
  
  if (!is.character(ylabel) && !is.na(ylabel))
     stop("'ylabel' input is incorrect, it should be of type character or string. Verify!")
  
  if (is.na(xlabel)) # || !is.character(xlabel))
     xlabel = "Effects"  # Nomeia Eixo X  
  
  if (is.na(ylabel)) # || !is.character(ylabel))
     ylabel = "Predicted values"  # Nomeia Eixo Y
  
  if (!is.logical(color))
     stop("'color' input is incorrect, it should be TRUE or FALSE. Verify!")
  
  if (!is.logical(casc))
     stop("'casc' input is incorrect, it should be TRUE or FALSE. Verify!")
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  Exp.Table <- table(BL$MPred[,1]) # tabela com as quantidade de amostras em cada experimento
  Exp.Names <- names(Exp.Table)    # nomes dos experimentos
  Num.Exp   <- length(Exp.Table)   # numero de experimentos

  if (Num.Exp != 0 && length(expcolor) != Num.Exp && !is.na(expcolor) ||
      Num.Exp == 0 && length(expcolor) != 1 && !is.na(expcolor))
    stop("'expcolor' input is incorrect, it should be in an amount equal to the number of experiments. Verify!")

  boxleg = ifelse(boxleg,"o","n") # moldura nas legendas, "n" sem moldura, "o" com moldura  
  
  cor <- 1 # cor inicial dos pontos e legendas
  ##### FIM - Informacoes usadas nos Graficos #####
  
  if (!is.character(titles[1]) || is.na(titles[1])) titles[1] = c("Study of the effects of concentrations")
  if (!is.character(titles[2]) || is.na(titles[2])) titles[2] = c("Component")
  
  Init.Form <- 15 # formato inicial dos pontos
  
  Form.Points <- Init.Form:(Init.Form + Num.Exp-1)
  
  if (color) {
    if (!is.na(expcolor[1])) {
      cor1 <- expcolor
    }
    else { cor1 <- cor:(cor + Num.Exp - 1) }
  }
  else { cor1 <- cor }
  
  #### INICIO - Efeitos da concentracao ####
  if (posleg == 1) Pos = "left" # posicao das legendas nos graficos
  if (posleg == 2) Pos = "right"
  if (posleg == 3) Pos = "top"
  if (posleg == 4) Pos = "bottom"
  
  Concentration <- BL$MPred[,3]
  
  if (casc) dev.new() # efeito cascata na apresentacao dos graficos
  
  plt <- xyplot(BL$MPred[,4] ~ BL$MPred[,2] | Concentration,
         xlab  = xlabel,
         ylab  = ylabel,
         main  = titles[1],
         data  = BL$MPred, # dados
         group = BL$MPred[,1], # grupos dos experimentos 
         type  = "b", # tipo de grafico
         col   = cor1,
         pch   = Form.Points,
         grid  = TRUE,
         as.table = TRUE, # desenha graficos da esquerda para direita de cima para baixo
         key = list(space = Pos, # criando a legenda
                    points = list(col = cor1, pch = Form.Points),
                    columns = ifelse(posleg < 3, 1, Num.Exp),
                    bty = boxleg,
                    text = list(paste("Exp -", Exp.Names))))
         # auto.key=list(space="right") # legenda automatica
  print(plt, split=c(1,1,1,1), more = F)
  #### FIM - Efeitos da concentracao ####

  
  #### INICIO - Preditos por componentes ####
  Data <- BL$MCPred
  Tit  <- colnames(Data[,3:ncol(Data)])

  # windows(width = 5, height = 5, pointsize = 1)

  nc <- (ncol(Data)-2) # numero de variaveis regressoras
  
  if (casc) dev.new() # efeito cascata na apresentacao dos graficos
  
  for(i in 1:nc) { # variaveis regressoras
    
    plt <- xyplot(Data[, 2 + i] ~ Data[,2],
           xlab  = xlabel,
           ylab  = ylabel,
           main  = paste(titles[2], Tit[i]),
           data  = Data, # dados
           group = Data[,1], # grupos dos experimentos
           type  = "b", # tipo de grafico
           col   = cor1,
           pch   = Form.Points,
           grid  = TRUE,
           as.table = TRUE, # desenha graficos da esquerda para direita de cima para baixo
           key = list(space = Pos, # criando a legenda
                      points = list(col = cor1, pch = Form.Points),
                      columns = ifelse(posleg < 3, 1, Num.Exp),
                      bty = boxleg,
                      text = list(paste("Exp -", Exp.Names))))
 
    num <- i/2
    if ((ceiling(num) - num) != 0) { # se for impar
        print(plt, split=c(1,1,1,2), more = ifelse(i < nc, TRUE, FALSE)) 
    } else {
        print(plt, split=c(1,2,1,2), more = FALSE)
        if (i < nc) if (casc) dev.new() # efeito cascata na apresentacao dos graficos
    }
  }
  #### FIM - Preditos por componentes ####
}
