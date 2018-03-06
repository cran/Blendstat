Plot.Blend <- function(BL, Titles = c(NA,NA), PosLeg = 2, xlabel = NA, 
                       ylabel = NA, BoxLeg = FALSE, Color = TRUE) {
  
  # Rotina para plotar graficos de Blendstate,
  # Desenvolvida por Marcelo Angelo Cirillo e
  # Paulo Cesar Ossani em 11/2017 
  
  # Entrada:
  # BL       - Dados da funcao Blend.
  # Titles   - Titulos para o grafico dos efeitos das concentracoes e componentes. Se nao for definido assume texto padrao.
  # PosLeg   - 1 para legenda a esquerda,
  #            2 para legenda a direita (default),
  #            3 para legenda acima,
  #            4 para legenda abaixo.
  # xlabel	 - Nomeia o eixo X, se nao definido retorna padrao.
  # ylabel	 - Nomeia o eixo Y, se nao definido retorna padrao.
  # BoxLeg   - Colocar moldura na legenda (default = TRUE).
  # Color    - Graficos coloridos (default = TRUE).

  # Retorna:
  # Varios graficos.
  
  if (!is.numeric(PosLeg) || PosLeg < 1 || PosLeg > 4 || (floor(PosLeg)-PosLeg) != 0)
     stop("Input to set the position of the legend 'PosLeg' is incorrect, must be an integer number between [1,4]. Check!")
  
  if (!is.logical(BoxLeg)) 
     stop("Input to insert the frame of the legend 'BoxLeg' is incorrect, must be TRUE or FALSE. Check!")
  
  if (!is.logical(Color))
     stop("Input for 'Color' is incorrect, must be TRUE or FALSE. Check!")

  if (is.na(xlabel) || !is.character(xlabel))
     xlabel = "Effects"  # Nomeia Eixo X  
  
  if (is.na(ylabel) || !is.character(ylabel))
     ylabel = "Predicted values"  # Nomeia Eixo Y
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  BoxLeg = ifelse(BoxLeg,"o","n") # moldura nas legendas, "n" sem moldura, "o" com moldura
  
  Exp.Table <- table(BL$MPred[,1]) # tabela com as quantidade de amostras em cada experimento
  Exp.Names <- names(Exp.Table)  # nomes dos experimentos
  Num.Exp   <- length(Exp.Table) # numero de experimentos

  cor <- 1 # cor inicial dos pontos e legendas
  ##### FIM - Informacoes usadas nos Graficos #####
  
  if (!is.character(Titles[1]) || is.na(Titles[1])) Titles[1] = c("Study of the effects of concentrations")
  if (!is.character(Titles[2]) || is.na(Titles[2])) Titles[2] = c("Component")
  
  Init.Form <- 15 # formato inicial dos pontos
  
  Form.Points <- Init.Form:(Init.Form + Num.Exp-1)
  
  if (Color) { 
     cor1 = cor:Num.Exp 
  } else cor1 = "black"
  
  #### INICIO - Efeitos da concentracao ####
  if (PosLeg == 1) Pos = "left" # posicao das legendas nos graficos
  if (PosLeg == 2) Pos = "right"
  if (PosLeg == 3) Pos = "top"
  if (PosLeg == 4) Pos = "bottom"
  
  Concentration <- BL$MPred[,3]
  
  dev.new() # nova tela para o grafico
  
  plt <- xyplot(BL$MPred[,4] ~ BL$MPred[,2] | Concentration,
         xlab  = xlabel,
         ylab  = ylabel,
         main  = Titles[1],
         data  = BL$MPred, # dados
         group = BL$MPred[,1], # grupos dos experimentos 
         type  = "b", # tipo de grafico
         col   = cor1,
         pch   = Form.Points,
         grid  = TRUE,
         # subscripts = "Concentracao",
         # layout = c(1,1),
         # jumper = F,
         # skip = F,
         # page = 1,
         as.table = TRUE, # desenha graficos da esquerda para direita de cima para baixo
         key = list(space = Pos, # criando a legenda
                    points = list(col = cor1, pch = Form.Points),
                    # title = "Experimentos",
                    columns = ifelse(PosLeg < 3, 1, Num.Exp),
                    bty = BoxLeg,
                    text = list(paste("Exp -", Exp.Names))))
         # auto.key=list(space="right") # legenda automatica
  print(plt, split=c(1,1,1,1), more = F)
  #### FIM - Efeitos da concentracao ####

  
  #### INICIO - Preditos por componentes ####
  Data <- BL$MCPred
  Tit  <- colnames(Data[,3:ncol(Data)])

  # windows(width = 5, height = 5, pointsize = 1)

  nc <- (ncol(Data)-2) # numero de variaveis regressoras
  
  dev.new() # nova tela para o grafico
  
  for(i in 1:nc) { # variaveis regressoras
    
    plt <- xyplot(Data[, 2 + i] ~ Data[,2],
           xlab  = xlabel,
           ylab  = ylabel,
           main  = paste(Titles[2], Tit[i]),
           data  = Data, # dados
           group = Data[,1], # grupos dos experimentos
           type  = "b", # tipo de grafico
           col   = cor1,
           pch   = Form.Points,
           grid = TRUE,
           as.table = TRUE, # desenha graficos da esquerda para direita de cima para baixo
           key = list(space = Pos, # criando a legenda
                      points = list(col = cor1, pch = Form.Points),
                      # title = "Experimentos",
                      columns = ifelse(PosLeg < 3, 1, Num.Exp),
                      bty = BoxLeg,
                      text = list(paste("Exp -", Exp.Names))))
 
    num <- i/2
    if ((ceiling(num) - num) != 0) { # se for impar
        print(plt, split=c(1,1,1,2), more = ifelse(i < nc, TRUE, FALSE)) 
    } else {
      print(plt, split=c(1,2,1,2), more = FALSE)
      if (i < nc) dev.new() # nova tela para o grafico
    }
  }
  #### FIM - Preditos por componentes ####
}
