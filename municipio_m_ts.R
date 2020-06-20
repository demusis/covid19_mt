# Automatizacao de processamentos relativos ao COVID-19
# Relatorio por municipio
# Autor: Carlo Ralph De Musis

# https://cran.r-project.org/web/views/TimeSeries.html
# https://www.datacamp.com/tracks/time-series-with-r
# https://www.datascience.com/blog/introduction-to-forecasting-with-arima-in-r-learn-data-science-tutorials
# https://www.otexts.org/fpp/8/7
# https://datascienceplus.com/time-series-analysis-using-arima-model-in-r/

# install.packages('robust')
# install.packages('TSPred')
# install.packages('hexbin')
# install.packages('tseries')
# install.packages('forecast')
# install.packages('fUnitRoots')
# install.packages('portes')
# install.packages('nortest')
# install.packages('tsoutliers')
# install.packages('boot')
# install.packages('imputeTS')
# install.packages('dygraphs')
# install.packages('changepoint')
# install.packages('drc')
# install.packages('envalysis')
# install.packages('sweep')
# install.packages('twitteR')
# install.packages('EpiEstim')
# install.packages('R0')
# install.packages('Metrics')
# install.packages('ggplot2')
# install.packages('dplyr')
# install.packages('hrbrthemes')
# install.packages('ggplot2')
# install.packages('reshape2')
# install.packages("EpiModel", dependencies = TRUE)
# install.packages("tictoc")
# install.packages('xlsx')
# install.packages('openxlsx')
# install.packages('deSolve')
# install.packages('network')
# install.packages('lubridate')
# install.packages('smooth')
# install.packages('Mcomp')
# install.packages('pracma')
# install.packages('units')
# install.packages('sf')

library(sf)
library(tidyverse)
library(RColorBrewer)
library(smooth)
library(Mcomp)
library(pracma)
#library(psych)
#library(dplyr)
#library(zoo)
suppressMessages(library(EpiModel))
library(tictoc)
library(lubridate)
# library(ggplot2)
# library(reshape2)
# library(xlsx)
library(openxlsx)
# library(Metrics)
library(EpiEstim)
# library(R0)
# library(TSPred)
# library(hexbin)
# library(robust)
library(imputeTS)
library(tseries)
library(forecast)
# library(fUnitRoots)
# library(portes)
# library(nortest)
# library(tsoutliers)
# library(boot)
# library(imputeTS)
# library(dygraphs)
# library(sweep)
# library(changepoint)
library(ggplot2)
library(drc)
# library(envalysis)
library(twitteR)
library(gridExtra)
library(zoo)
library(dplyr)
library(hrbrthemes)
library(Metrics)

# Parametros gerais VA.
peso_infectados_va <- 1
peso_enfermarias_va <- 0
peso_uti_va <- 0

peso_enfermarias_to <- 0
peso_uti_to <- 1

# Liga o cronometro
tic("Inicio...")

df_lolipop <- data.frame(
  # ers = character(), ## Nome do Escritorio Regional de Saude.
  municipio = character(), ## Nome do Municipio/Estado.
  infectados_n = character(), ## Numero de infectados.
  recuperados_n  = character(), ## Numero de recuperados.
  obitos_n = character(), ## Numero de obitos.
  r_025 = double(), ## Percentil de 2.5% para o Rt.
  r_500 = double(), ## Mediana do Rt.
  r_975 = double(), ##  Percentil de 97.5% para o Rt.
  r_mm7  = double(), # Media movel de 7 dias de r_500
  n_va = double(), ## Velocidade de avanco de casos novos.
  p_va = double(), ## Velocidade de avanco ponderada.
  to_uti = double(), # Taxa de ocupação de UTI.
  stringsAsFactors = FALSE
)

df_estado <- data.frame(
  r_025 = double(), # Percentil de 2.5% para o Rt.
  r_500 = double(), # Mediana do Rt.
  r_975 = double(), #  Percentil de 97.5% para o Rt. 
  p_va = double(), # Velocidade de avanco ponderada.
  n_va = double(), ## Velocidade de avanco de casos novos.
  enfermarias_va = double(), # Velocidade de avanco de enfermarias.
  uti_va = double(), # Velocidade de avanco de UTIs.
  to_uti = double(), # Taxa de ocupação de UTI.  
  stringsAsFactors = FALSE
)

# Define e cria dataframe para parametros epidemiologicos por municipio. 
par_df_sim <- data.frame(
  municipio = character(),
  inf.prob = double(),
  act.rate = double(),
  dias_rec = double(),
  mse = double(),
  stringsAsFactors = FALSE
)

# Define e cria dataframe para projecoes por aglomerado/municipio.
proj_aglomerado <- data.frame(
  Aglomerado = character(),
  'Dia Juliano' = double(),
  Previsao  = double(),
  'Limite superior (95%)'  = double(),
  stringsAsFactors = FALSE  
)

# Configuracao do Twitter
#try(setup_twitter_oauth(
#    consumer_key = "K8H2usTrUWKWrdUDhqZoMtYxT",
#    access_token = "1254974121265500160-ZW0Gy3GLNe1xWsMEasMUzQ0G67rufK",
#    consumer_secret = "egNsRyTB1qLJciGT6qweJjFhc2ZOXdU4t86Zfx9x0G9hD1rgUa",
#    access_secret = "YKNvlxrp5VhgVhJ3zZpP5KFalZPkKzC4VtnKS56JUoxtb"
#))

# Define diretorio de trabalho e carrega progresso de infectados
setwd("D:/Mega/COVID-19")
# setwd("~/COVID-19")

# Carrega aglomerados
# aglomerados <- read.xlsx("aglomerados.xlsx", 1, header = TRUE)
aglomerados <-
  read.csv("aglomerados.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

# Apaga eventuais celulas vazias.
aglomerados <- aglomerados[!is.na(aglomerados$aglomerado), ]

# Converte para data
aglomerados$data <-
  as.Date(parse_date_time(aglomerados[["data"]], '%d%m%y'), format = '%d%m%y')

# Para resolver eventuais problemas de conversao
# aglomerados[2:6] <- lapply(aglomerados[2:6], as.numeric)

# Altera nome de coluna, para padronizacao.
colnames(aglomerados)[1] <- 'aglomerado'

# Abre arquivo de parametros para o Mato Grosso.
# 80-90% nÃ£o precisam de internacao.
# 20% dos internados em UTI morrem.
# 5-15% dos internados precisam de UTI.
total_param_mt <-
  read.csv("mato_grosso.CSV", stringsAsFactors = FALSE)
total_param_mt$data <- dmy(total_param_mt$data)
total_param_mt$Dia_Juliano <- yday(total_param_mt$data)

# Carrega as series temporais da epidemia.
aux_dados <- read.csv("evolucao_ers.CSV", stringsAsFactors = FALSE)
aux_dados$Data <- dmy(aux_dados$Data)
aux_dados$Dia_Juliano <- yday(aux_dados$Data)
aux_dados$Total <-
  aux_dados$Infectados + aux_dados$Recuperados + aux_dados$Obitos

# Apaga eventuais celulas vazias.
aux_dados <- aux_dados[!is.na(aux_dados[, 'Infectados']),]

# Seleciona campos.
mt_dados <- aux_dados[, c('Municipio', 'Dia_Juliano', 
                          'Infectados', 'Recuperados',
                          'Obitos', 'Total')]

# Agrega os dados por ERS. # Desnecessario?
# mt_aux_dados <- aggregate(Total ~ Municipio + Dia_Juliano, mt_dados, sum)
mt_aux_dados <- mt_dados[, c('Total', 'Municipio', 'Dia_Juliano')]

# Agrega todos os dados para o Mato Grosso.
mt_total_dados <- aggregate(Total ~ Dia_Juliano, mt_dados, sum)
mt_total_dados$Municipio <- 'Mato Grosso'

# Adiciona o Mato Grosso no rol.
mt_aux_dados <- rbind(mt_aux_dados, mt_total_dados)

# Cria fatores
mt_aux_dados$Municipio <- factor(mt_aux_dados$Municipio)


#
# Funcao velocidade de avanco em 1 semana
#
va <- function(or_st) {
  st <- data.frame(dj = integer(),
                   y = double(),
                   stringsAsFactors = FALSE)
  posicao <- end(or_st)[1] - start(or_st)[1] + 1
  for (n_posicao in 15:posicao) {
    valor_t <- or_st[n_posicao]
    valor_t_1 <- or_st[n_posicao - 7][1]
    valor_t_2 <- or_st[n_posicao - 14][1]
    aux_st <- data.frame(time(or_st)[n_posicao],
                         (valor_t - valor_t_1) / (valor_t_1 - valor_t_2) - 1)
    colnames(aux_st) <- c('dj', 'y')
    st <- rbind(st, aux_st)
  }
  aux_st <- ts(st$y,
               start = st[1, 'dj'],
               frequency = 1)
  aux_st[is.nan(aux_st)] <- NA
  aux_st[is.infinite(aux_st)] <- NA
  return(aux_st)
}


# Resgata data do ultimo registro.
aux_dh <- aux_dados[nrow(aux_dados), 'Data']

# Laco para todos os municipios e Mato Grosso
for (aux_nivel in levels(mt_aux_dados$Municipio)) {
  # Seleciona dados do municipio.
  municipio_dados <- mt_aux_dados[mt_aux_dados$Municipio == aux_nivel,]
  seq_dados <-
    data.frame(seq(
      min(municipio_dados$Dia_Juliano),
      max(municipio_dados$Dia_Juliano),
      by = 1
    ))
  colnames(seq_dados) <- c('Dia_Juliano')
  
  # Cria serie numerica fechada para verificar ausentes. 
  todos_dados <-
    merge(x = seq_dados,
          y = municipio_dados,
          by = 'Dia_Juliano',
          all.x = TRUE)
  todos_dados$Municipio <- aux_nivel
  
  # Seleciona variaveis de analise.
  dados <- todos_dados[, c('Dia_Juliano', 'Total')]
  
  # Cria objeto TS.
  infectados <-
    ts(dados$Total,
       start = dados[1, 'Dia_Juliano'],
       frequency = 1)
  
  # Preenchimento de falhas
  infectados <-
    na_kalman(infectados, model = "StructTS", smooth = TRUE)
  
  # Insere dados iniciais na tabela de municipios.
  if (aux_nivel!="Mato Grosso") {
    aux_df_lolipop <- df_lolipop[0, ]
    
    mun_aux_dados <- aux_dados[aux_dados$Municipio == aux_nivel,]
    mun_aux_dados <-
      mun_aux_dados[which.max(mun_aux_dados$Dia_Juliano),]
    
    # aux_df_lolipop[1, 'ers'] <- mun_aux_dados$ERS
    aux_df_lolipop[1, 'municipio'] <- mun_aux_dados$Municipio
    aux_df_lolipop[1, 'infectados_n'] <- mun_aux_dados$Infectados
    aux_df_lolipop[1, 'recuperados_n'] <- mun_aux_dados$Recuperados
    aux_df_lolipop[1, 'obitos_n'] <- mun_aux_dados$Obitos
  }
  # Filtra as series com ao menos 15 dias de dados e 10 infectados.
  if (max(infectados) >= 10 && length(infectados) >= 15) {
    print(paste('TS - municipio -', aux_nivel))
    
    # Pequisa transformacao de Box-Cox
    l <- BoxCox.lambda(infectados)
    
    # Testa os modelos disponiveis no pacote Forecast
    try({
      ets_ajuste <- ets(infectados, lambda = l)
      arima_ajuste <-
        auto.arima(
          infectados,
          stepwise = FALSE,
          approximation = FALSE,
          trace = FALSE,
          lambda = l
        )
      bats_ajuste <- bats(infectados, lambda = l)
      bagged_ajuste <- baggedModel(infectados, lambda = l)
      croston_ajuste <- croston(infectados)
      tslm_ajuste <- tslm(infectados ~ trend)
      # try(nnetar_ajuste <- nnetar(infectados, repeats = 10, lambda = l))
      tbats_ajuste <- tbats(infectados, use.box.cox = TRUE)
      theta_ajuste <- thetaf(infectados)
    })
    
    # Verifica qual dos modelos foi o mais acurado
    modelo <- ets_ajuste
    acuracia_modelo <- forecast::accuracy(modelo)[2]
    for (aux_teste in list(
      arima_ajuste,
      bats_ajuste,
      bagged_ajuste,
      croston_ajuste,
      theta_ajuste,
      tslm_ajuste,
      # nnetar_ajuste, # O modelo por RN apresenta alguns problemas para validacao
      tbats_ajuste
    )) {
      try({
        aux_acuracia <- forecast::accuracy(aux_teste)[2]
        if (aux_acuracia < acuracia_modelo) {
          modelo <- aux_teste
          acuracia_modelo <- aux_acuracia
        }
      })
    }
    
    # Previso para 15 dias (solicitacao da SES).
    try(previsao <- forecast(modelo, level = c(90), h = 15))
    df_previsao <- data.frame(previsao)
    colnames(df_previsao) <-
      c('Previsao', 'Limite inferior (95%)', 'Limite superior (95%)')
    
    # Taxas de crescimento.
    df_previsao$Aglomerado <- aux_nivel
    df_previsao <- tibble::rownames_to_column(df_previsao, var = "Dia Juliano")
    proj_aglomerado <- rbind(proj_aglomerado, df_previsao)
    
    # Velocidade de avanco de infectados.
    ts_in1s <- va(infectados)
    
    if ((end(ts_in1s)[1] - start(ts_in1s)[1]) >= 14) {
      try(modelo_in1s <- sma(
        tail(ts_in1s, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = 'all'
      ))
      
      n_va <- last(modelo_in1s$forecast)
    } else {
      n_va <- NA
    }
    
    # Estima R
    dados_MT <- as.numeric(diff(infectados))
    if (min(dados_MT)>=0) {
    res <- NULL
    try(res <- estimate_R(dados_MT, method = "parametric_si",
                          config = make_config(list(
                            mean_si = 3.96, std_si = 4.75
                          ))))
    if (!is.null(res)) {
      r_025 <- tail(res$R$'Quantile.0.025(R)', 1)
      r_975 <- tail(res$R$'Quantile.0.975(R)', 1)
      r_500 <- tail(res$R$'Median(R)', 1)
      
      # Avalia a media movel do R dos ultimos 30 dias.
      r30 <- tail(res$R$'Median(R)', 30)
      mm_r30 <- sma(
        r30,
        order = 7,
        h = 1,
        interval = 'none',
        silent = "all"
      )
      Rt <- mm_r30$forecast
    } else {
      r_025 <- NA
      r_975 <- NA
      r_500 <- NA
      Rt <- NA
    }
    } else {
      r_025 <- -999
      r_975 <- -999
      r_500 <- -999
      Rt <- -999
      
    }
    
    if (aux_nivel != 'Mato Grosso') {
      # Insere estatísticas dos municipios.
      aux_df_lolipop[1, 'r_025'] <- r_025
      aux_df_lolipop[1, 'r_500'] <- r_500
      aux_df_lolipop[1, 'r_975'] <- r_975
      aux_df_lolipop[1, 'r_mm7'] <- Rt
      aux_df_lolipop[1, 'n_va'] <- n_va
      aux_df_lolipop[1, 'p_va'] <- NA # Depende dos dados de MT.
    }
    # Estima estatistica do Estado
    else
    {
      # Velocidade de avanco em enfermarias utilizadas
      enfermarias <-
        ts(
          total_param_mt$enfermaria_utilizadas,
          start = total_param_mt[1, 'Dia_Juliano'] - 1,
          frequency = 1
        )
      
      ts_enfermarias <- va(enfermarias)
      modelo_enfermarias <- sma(
        tail(ts_enfermarias, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = 'all'
      )

      # Velocidade de avanco em UTIs
      UTI <-
        ts(total_param_mt$uti_utilizadas,
           start = total_param_mt[1, 'Dia_Juliano'] - 1,
           frequency = 1)
      
      ts_UTI <- va(UTI)
      modelo_UTI <- sma(
        tail(ts_UTI, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = 'all'
      )

      # Velocidade de velocidade de avanco ponderada.
      pva_geral <- (as.numeric(modelo_in1s$y)*peso_infectados_va +
                    as.numeric(modelo_enfermarias$y)*peso_enfermarias_va +
                    as.numeric(modelo_UTI$y)*peso_uti_va)/
                    (peso_infectados_va + peso_enfermarias_va + peso_uti_va) 
      ts_pva_geral <-
        ts(pva_geral,
           start = total_param_mt[1, 'Dia_Juliano'] - 1,
           frequency = 1)
      modelo_pva <- sma(
        tail(ts_pva_geral, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = 'all'
      )
      df_estado[1, 'r_025'] <- r_025
      df_estado[1, 'r_500'] <- r_500
      df_estado[1, 'r_975'] <- r_975
      df_estado[1, 'r_mm7'] <- Rt
      df_estado[1, 'p_va'] <- modelo_pva$forecast  
      df_estado[1, 'n_va'] <- n_va
      df_estado[1, 'enfermarias_va'] <- modelo_enfermarias$forecast
      df_estado[1, 'uti_va'] <- modelo_UTI$forecast  
      
      df_estado[1, 'uti_to'] <- total_param_mt[which.max(total_param_mt$Dia_Juliano), 'to_leitos_uti_sus']
      df_estado[1, 'enfermarias_to'] <- total_param_mt[which.max(total_param_mt$Dia_Juliano), 'to_leitos_clinicos_sus']
      df_estado[1, 'iv'] <- (df_estado[1, 'enfermarias_to']*peso_enfermarias_to + 
                                df_estado[1, 'uti_to']*peso_uti_to)/
                                (peso_enfermarias_va + peso_uti_to)
      df_estado[1, 'infectados_n'] <- total_param_mt[which.max(total_param_mt$Dia_Juliano), 'casos_confirmados']
      df_estado[1, 'recuperados_n']  <- total_param_mt[which.max(total_param_mt$Dia_Juliano), 'casos_recuperados']
      df_estado[1, 'obitos_n'] <- total_param_mt[which.max(total_param_mt$Dia_Juliano), 'obitos']
      
      df_estado$iv_cat <- NA
      if (df_estado[1, 'iv']<.30) {
        df_estado[1, 'iv_cat'] <- "Muito baixo"
      } else {
        if (df_estado[1, 'iv']<.60) {
          df_estado[1, 'iv_cat'] <- "Baixo"
        } else {
          if (df_estado[1, 'iv']<.75) {
            df_estado[1, 'iv_cat'] <- "Moderado"
          } else {
            if (df_estado[1, 'iv']<.90) {
              df_estado[1, 'iv_cat'] <- "Alto"
            } else {
              df_estado[1, 'iv_cat'] <- "Muito alto"
            }
          }
        }
      }
    }
    rm(r_025, r_500, r_975, Rt, 
       modelo_in1s, modelo_pva, modelo_enfermarias, modelo_UTI,
       n_va, res
       )
  }
  if (aux_nivel!="Mato Grosso") {
    df_lolipop <- rbind(df_lolipop, aux_df_lolipop)
  }
}

df_lolipop$municipio <- factor(df_lolipop$municipio)
for (aux_nivel in levels(df_lolipop$municipio)) {
  municipio_df_lolipop <- df_lolipop[df_lolipop$municipio==aux_nivel, ]
  if (!is.na(df_lolipop[df_lolipop$municipio==aux_nivel, 'n_va'])) {
    df_lolipop[df_lolipop$municipio==aux_nivel, 'p_va'] <-
      (municipio_df_lolipop$n_va*peso_infectados_va + 
       df_estado$enfermarias_va*peso_enfermarias_va + 
       df_estado$uti_va*peso_uti_va)/
      (peso_infectados_va + peso_enfermarias_va + peso_uti_va)
  } else{
    df_lolipop[df_lolipop$municipio==aux_nivel, 'p_va'] <- NA
  }
}

ibge_municipios <-
  read.csv("ibge_municipios.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

df_lolipop <- merge(ibge_municipios, df_lolipop, all.x = TRUE)
df_lolipop[is.na(df_lolipop$infectados_n), 'infectados_n'] <- 0
df_lolipop[is.na(df_lolipop$recuperados_n), 'recuperados_n'] <- 0
df_lolipop[is.na(df_lolipop$obitos_n), 'obitos_n'] <- 0

setwd("D:/Mega/COVID-19/resultados")


#
# Classifica conforme as ameacas
#

cIA <- function(transmissao="TRANSMISSAO COMUNITARIA", infectados, Rt, Va) {
  aux_risco <- NA
  if (is.na(Rt)) {
    Rt <- -999
  }
  if (is.na(Va)) {
    Va <- -999
  }
  if (infectados == 0 ||
      ((transmissao == "TRANSMISSAO LOCAL" ||
        transmissao == "CASOS ISOLADOS") &&
      Rt < 1 &&
      Va < 0.15)) {
    aux_risco <- "Muito baixo"
  }
  if ((transmissao == "TRANSMISSAO LOCAL" ||
       transmissao == "CASOS ISOLADOS") &&
      (Rt >= 1.15 ||
       Va >= 0.25)) {
    aux_risco <- "Baixo"
  }
  if (transmissao == "TRANSMISSAO COMUNITARIA" ||
      (Rt >= 1.30 &&
       Va >= 0.40)) {
    aux_risco <- "Moderado"
  }
  if (transmissao == "TRANSMISSAO COMUNITARIA" &&
      (Rt >= 1.60 &&
       Va >= 1)) {
    aux_risco <- "Alto"
  }
  if (transmissao == "TRANSMISSAO COMUNITARIA" &&
      (Rt >= 2 &&
       Va >= 1.7)) {
    aux_risco <- "Muito alto"
  }
  return(aux_risco)
}

df_lolipop$municipio <- factor(df_lolipop$municipio)
df_lolipop$ia_cat <- NA
for (aux_nivel in levels(df_lolipop$municipio)) {
  df_lolipop[df_lolipop$municipio == aux_nivel, 'ia_cat'] <-
    cIA(df_lolipop[df_lolipop$municipio == aux_nivel, 'classificacao_epidemiologica'],
        df_lolipop[df_lolipop$municipio == aux_nivel, 'infectados_n'],
        df_lolipop[df_lolipop$municipio == aux_nivel, 'r_500'],
        df_lolipop[df_lolipop$municipio == aux_nivel, 'p_va'])
}

df_estado$ia_cat <- NA
df_estado[1, 'ia_cat'] <-
  cIA("TRANSMISSAO COMUNITARIA",
      Inf,
      df_estado[1, 'r_mm7'],
      df_estado[1, 'p_va'])


# Matriz de risco
matriz_risco <- data.frame(
  ia_cat = character(),
  iv_cat = character(),
  risco = character(),
  stringsAsFactors = FALSE
)
# c1
matriz_risco[1, 'ia_cat'] <- "Muito alto" 
matriz_risco[1, 'iv_cat'] <- "Muito baixo"
matriz_risco[1, 'risco'] <- "Muito alto"

matriz_risco[2, 'ia_cat'] <- "Alto" 
matriz_risco[2, 'iv_cat'] <- "Muito baixo"
matriz_risco[2, 'risco'] <- "Alto"

matriz_risco[3, 'ia_cat'] <- "Moderado" 
matriz_risco[3, 'iv_cat'] <- "Muito baixo"
matriz_risco[3, 'risco'] <- "Moderado"

matriz_risco[4, 'ia_cat'] <- "Baixo" 
matriz_risco[4, 'iv_cat'] <- "Muito baixo"
matriz_risco[4, 'risco'] <- "Baixo"

matriz_risco[5, 'ia_cat'] <- "Muito baixo" 
matriz_risco[5, 'iv_cat'] <- "Muito baixo"
matriz_risco[5, 'risco'] <- "Muito baixo"

# c2
matriz_risco[6, 'ia_cat'] <- "Muito alto" 
matriz_risco[6, 'iv_cat'] <- "Baixo"
matriz_risco[6, 'risco'] <- "Muito alto"

matriz_risco[7, 'ia_cat'] <- "Alto" 
matriz_risco[7, 'iv_cat'] <- "Baixo"
matriz_risco[7, 'risco'] <- "Alto"

matriz_risco[8, 'ia_cat'] <- "Moderado" 
matriz_risco[8, 'iv_cat'] <- "Baixo"
matriz_risco[8, 'risco'] <- "Moderado"

matriz_risco[9, 'ia_cat'] <- "Baixo" 
matriz_risco[9, 'iv_cat'] <- "Baixo"
matriz_risco[9, 'risco'] <- "Moderado"

matriz_risco[10, 'ia_cat'] <- "Muito baixo" 
matriz_risco[10, 'iv_cat'] <- "Baixo"
matriz_risco[10, 'risco'] <- "Baixo"

# c3
matriz_risco[11, 'ia_cat'] <- "Muito alto" 
matriz_risco[11, 'iv_cat'] <- "Moderado"
matriz_risco[11, 'risco'] <- "Muito alto"

matriz_risco[12, 'ia_cat'] <- "Alto" 
matriz_risco[12, 'iv_cat'] <- "Moderado"
matriz_risco[12, 'risco'] <- "Alto"

matriz_risco[13, 'ia_cat'] <- "Moderado" 
matriz_risco[13, 'iv_cat'] <- "Moderado"
matriz_risco[13, 'risco'] <- "Alto"

matriz_risco[14, 'ia_cat'] <- "Baixo" 
matriz_risco[14, 'iv_cat'] <- "Moderado"
matriz_risco[14, 'risco'] <- "Moderado"

matriz_risco[15, 'ia_cat'] <- "Muito baixo" 
matriz_risco[15, 'iv_cat'] <- "Moderado"
matriz_risco[15, 'risco'] <- "Moderado"

# c4
matriz_risco[16, 'ia_cat'] <- "Muito alto" 
matriz_risco[16, 'iv_cat'] <- "Alto"
matriz_risco[16, 'risco'] <- "Muito alto"

matriz_risco[17, 'ia_cat'] <- "Alto" 
matriz_risco[17, 'iv_cat'] <- "Alto"
matriz_risco[17, 'risco'] <- "Muito alto"

matriz_risco[18, 'ia_cat'] <- "Moderado" 
matriz_risco[18, 'iv_cat'] <- "Alto"
matriz_risco[18, 'risco'] <- "Alto"

matriz_risco[19, 'ia_cat'] <- "Baixo" 
matriz_risco[19, 'iv_cat'] <- "Alto"
matriz_risco[19, 'risco'] <- "Alto"

matriz_risco[20, 'ia_cat'] <- "Muito baixo" 
matriz_risco[20, 'iv_cat'] <- "Alto"
matriz_risco[20, 'risco'] <- "Alto"

# c5
matriz_risco[21, 'ia_cat'] <- "Muito alto" 
matriz_risco[21, 'iv_cat'] <- "Muito alto"
matriz_risco[21, 'risco'] <- "Muito alto"

matriz_risco[22, 'ia_cat'] <- "Alto" 
matriz_risco[22, 'iv_cat'] <- "Muito alto"
matriz_risco[22, 'risco'] <- "Muito alto"

matriz_risco[23, 'ia_cat'] <- "Moderado" 
matriz_risco[23, 'iv_cat'] <- "Muito baixo"
matriz_risco[23, 'risco'] <- "Muito alto"

matriz_risco[24, 'ia_cat'] <- "Baixo" 
matriz_risco[24, 'iv_cat'] <- "Muito alto"
matriz_risco[24, 'risco'] <- "Muito alto"

matriz_risco[25, 'ia_cat'] <- "Muito baixo" 
matriz_risco[25, 'iv_cat'] <- "Muito alto"
matriz_risco[25, 'risco'] <- "Muito alto"

matriz_risco$ia_cat <- factor(matriz_risco$ia_cat)
matriz_risco$iv_cat <- factor(matriz_risco$iv_cat)

# Utiliza a categoria de vulnerabilidade do Estado nos municípios.
df_lolipop$iv_cat <- df_estado[1, 'iv_cat']
df_lolipop$to_uti <- df_estado[1, 'iv']

# Classifica os riscos
df_lolipop$risco <- NA
df_estado$risco <- NA
for (aux_ia in levels(matriz_risco$ia_cat)) {
  for (aux_iv in levels(matriz_risco$iv_cat)) {
    aux_risco <- matriz_risco[((matriz_risco$ia_cat==aux_ia) &
                               (matriz_risco$iv_cat==aux_iv)), 'risco']
    df_lolipop[((df_lolipop$ia_cat == aux_ia) & 
                (df_lolipop$iv_cat == aux_iv) &
                !is.na(df_lolipop$ia_cat) &
                !is.na(df_lolipop$iv_cat)), 'ia_cat'] -> ia_aux_valor

    if (!identical(ia_aux_valor, character(0))) {
      df_lolipop[((df_lolipop$ia_cat == aux_ia) &
                 (df_lolipop$iv_cat == aux_iv) &
                  !is.na(df_lolipop$ia_cat) &
                  !is.na(df_lolipop$iv_cat)), 'risco'] <- aux_risco
    }
    df_estado[df_estado$ia_cat == aux_ia, 'ia_cat'] -> ia_aux_valor
    df_estado[((df_estado$ia_cat == aux_ia) & 
               (df_estado$iv_cat == aux_iv) &
                !is.na(df_estado$ia_cat) &
                !is.na(df_estado$iv_cat)), 'ia_cat'] -> ia_aux_valor  
    if (!identical(ia_aux_valor, character(0))) {
      df_estado[((df_estado$ia_cat == aux_ia) & 
                 (df_estado$iv_cat == aux_iv) &
                  !is.na(df_estado$ia_cat) &
                  !is.na(df_estado$iv_cat)), 'risco'] <- aux_risco
    }
  }
}

# Salva as estatisticas por aglomerado
write.table(
  file = paste(aux_dh, '- bi_aglomerados.csv'),
  df_lolipop,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ';',
  dec = ',',
  na = ''
)

write.table(
  file = paste(aux_dh, '- bi_mt.csv'),
  df_estado,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ';',
  dec = ',',
  na = ''
)

#
# Roda classificacao de risco Fuzzy em Python
#

# install.packages('reticulate')
library(reticulate)

# Para verificar path do Python.
# >>> import os
# >>> import sys
# >>> os.path.dirname(sys.executable)

use_python('C:/Users/Carlo/AppData/Local/Programs/Python/Python37')
py_run_file('risco_fuzzy.py')

head(py$df)

