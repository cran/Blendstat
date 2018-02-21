Plot.Blend <- function(BL, Title = NULL, PosLeg = 2, BoxLeg = FALSE, Color = TRUE) {
  
  # Rotina para plotar graficos de Blendstate,
  # Desenvolvida por Marcelo Angelo Cirillo e
  # Paulo Cesar Ossani em 11/2017 
  
  # Entrada:
  # BL       - Dados da funcao Blend.
  # Title    - Titulos para o grafico dos efeitos das concentracoes. Se nao for definido assume texto padrao.
  # PosLeg   - 1 para legenda a esquerda,
  #            2 para legenda a direita (default),
  #            3 para legenda acima,
  #            4 para legenda abaixo.
  # BoxLeg   - Colocar moldura na legenda (default = TRUE).
  # Color    - Graficos coloridos (default = TRUE).
  
  # BL=Res; Title = NULL; PosLeg = 2; BoxLeg = FALSE; Color = TRUE
  
  # Retorna:
  # Varios graficos.
  
  if (!is.numeric(PosLeg) || PosLeg < 1 || PosLeg > 4 || (floor(PosLeg)-PosLeg) != 0)
    stop("Input to set the position of the legend 'PosLeg' is incorrect, must be an integer number between [1,4]. Check!")
  
  if (!is.logical(BoxLeg)) 
    stop("Input to insert the frame of the legend 'BoxLeg' is incorrect, must be TRUE or FALSE. Check!")
  
  if (!is.logical(Color))
    stop("Input for 'Color' is incorrect, must be TRUE or FALSE. Check!")

  ##### INICIO - Informacoes usadas nos Graficos #####
  BoxLeg = ifelse(BoxLeg,"o","n") # moldura nas legendas, "n" sem moldura, "o" com moldura
  
  Exp.Table <- table(BL$MPred[,1]) # tabela com as quantidade de amostras em cada experimento
  Exp.Names <- names(Exp.Table)  # nomes dos experimentos
  Num.Exp   <- length(Exp.Table) # numero de experimentos

  cor <- 1 # cor inicial dos pontos e legendas
  ##### FIM - Informacoes usadas nos Graficos #####
  
  if (!is.character(Title[1]) || is.na(Title[1])) Title[1] = c("Study of the effects of concentrations")
 
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
  
  plt <- xyplot(BL$MPred[,4] ~ BL$MPred[,2] | Concentration,
         xlab  = "Effects",
         ylab  = "Predicted values",
         main  = Title[1],
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

  for(i in 1:nc) { # variaveis regressoras
 
    plt <- xyplot(Data[, 2 + i] ~ Data[,2],
           xlab  = "Effects",
           ylab  = "Predicted values",
           main  = paste("Component", Tit[i]),
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
    } else print(plt, split=c(1,2,1,2), more = FALSE)

  }
  #### FIM - Preditos por componentes ####
}