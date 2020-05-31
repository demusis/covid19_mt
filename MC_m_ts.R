# Automatizacao de processamentos relativos ao COVID-19
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

# Carrega o EpiModel e demais bibliotecas
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
#library(envalysis)
library(twitteR)
library(gridExtra)
library(zoo)
library(dplyr)
library(hrbrthemes)
library(Metrics)

# Liga o cronometro
tic("Inicio...")

# Define e cria dataframe para grafico de lolipop.
df_lolipop <- data.frame(
  ers = character(),
  p025 = double(),
  p50 = double(),
  p975 = double(),
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

# Define e cria dataframe para projecoe por aglomerado/municipio.
proj_aglomerado <- data.frame(
  Aglomerado = character(),
  'Dia Juliano' = double(),
  Previsao  = double(),
  'Limite superior (95%)'  = double(),
  stringsAsFactors = FALSE  
)

# Configuracao do Twitter
try(setup_twitter_oauth(
  consumer_key = "K8H2usTrUWKWrdUDhqZoMtYxT",
  access_token = "1254974121265500160-ZW0Gy3GLNe1xWsMEasMUzQ0G67rufK",
  consumer_secret = "egNsRyTB1qLJciGT6qweJjFhc2ZOXdU4t86Zfx9x0G9hD1rgUa",
  access_secret = "YKNvlxrp5VhgVhJ3zZpP5KFalZPkKzC4VtnKS56JUoxtb"
))

# Define diretorio de trabalho e carrega progresso de infectados
# setwd("D:/OneDrive/Notebooks/Python/Corona virus/SARIMA")
setwd("~/COVID-19")

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
# 80-90% não precisam de internacao.
# 20% dos internados em UTI morrem.
# 5-15% dos internados precisam de UTI.
total_param_mt <-
  read.csv("mato_grosso.CSV", stringsAsFactors = FALSE)
total_param_mt$data <- dmy(total_param_mt$data)
total_param_mt$Dia_Juliano <- yday(total_param_mt$data)

# Carrega as series temporais da epidemia para o Mato Grosso.
aux_dados <- read.csv("evolucao_ers.CSV", stringsAsFactors = FALSE)
aux_dados$Data <- dmy(aux_dados$Data)
aux_dados$Dia_Juliano <- yday(aux_dados$Data)
aux_dados$Total <-
  aux_dados$Infectados + aux_dados$Recuperados + aux_dados$Obitos

# Apaga eventuais celulas vazias.
aux_dados <- aux_dados[!is.na(aux_dados[, 'Infectados']),]

# Calcula e insere registro para o Mato Grosso.
aux_dados$ERS <- factor(aux_dados$ERS)
mt_dados <- aux_dados[, c('ERS', 'Dia_Juliano', 'Total')]
mt_aux_dados <- aggregate(Total ~ ERS + Dia_Juliano, mt_dados, sum)

mt_total_dados <- aggregate(Total ~ Dia_Juliano, mt_dados, sum)
mt_total_dados$ERS <- 'Mato Grosso'
mt_aux_dados <- rbind(mt_aux_dados, mt_total_dados)

# mt_aux_dados <- mt_aux_dados[, -2]
# colnames(mt_aux_dados) <- c('Dia_Juliano', 'Total')
# sum(aux_dados[,])

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
# aux_dh <- Sys.time()
aux_dh <- aux_dados[nrow(aux_dados), 'Data']

# Laco para todas ERS e Mato Grosso
for (aux_nivel in levels(mt_aux_dados$ERS)) {
  ers_dados <- mt_aux_dados[mt_aux_dados$ERS == aux_nivel,]
  seq_dados <-
    data.frame(seq(
      min(ers_dados$Dia_Juliano),
      max(ers_dados$Dia_Juliano),
      by = 1
    ))
  colnames(seq_dados) <- c('Dia_Juliano')
  
  todos_dados <-
    merge(x = seq_dados,
          y = ers_dados,
          by = 'Dia_Juliano',
          all.x = TRUE)
  todos_dados$ERS <- aux_nivel
  
  dados <- todos_dados[, c('Dia_Juliano', 'Total')]
  
  # Cria objeto TS
  infectados <-
    ts(dados$Total,
       start = dados[1, 'Dia_Juliano'],
       frequency = 1)
  
  # Preenchimento de falhas
  infectados <-
    na_kalman(infectados, model = "StructTS", smooth = TRUE)
  
  # Filtra as series com ao menos 15 dias de dados e 10 infectados
  if (max(infectados) >= 10 && length(infectados) >= 15) {
    print(paste('TS -', aux_nivel))
    
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
    
    # Previsao para 3 dias
    previsao <- forecast(modelo, level = c(95), h = 3)
    
    # Posta previsao no twitter
    jpeg(
      "grafico_temp.jpg",
      width = 1400,
      height = 700,
      quality = 90
    )
    
    plot(previsao,
         main = aux_nivel,
         xlab = "Dias",
         ylab = "Infectados")
    dev.off()
    
    # Posta tweets da previsao de curto prazo
    try(updateStatus(
      paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Previsões de curto prazo  [Gráfico - 3 dias].'
      ),
      mediaPath = "grafico_temp.jpg"
    ))
    
    # Previso para 15 dias (solicitação daa SES),
    try(previsao <- forecast(modelo, level = c(90), h = 15))
    jpeg(
      "tabela_temp.jpg",
      width = 400,
      height = 450,
      quality = 90
    )
    df_previsao <-
      data.frame(previsao)[, c('Point.Forecast', 'Hi.100')]
    colnames(df_previsao) <-
      c('Previsao', 'Limite superior (95%)')
    grid.table(df_previsao)
    dev.off()
    
    # Posta tabela no Twitter.
    try(updateStatus(
      paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Previsões de curto prazo [Tabela - 15 dias]'
      ),
      mediaPath = "tabela_temp.jpg"
    ))
    
    # Taxas de crescimento
    df_previsao$Aglomerado <- aux_nivel
    df_previsao <- tibble::rownames_to_column(df_previsao, var = "Dia Juliano")
    proj_aglomerado <- rbind(proj_aglomerado, df_previsao)
    
    # Taxa de crescimento diario
    dif_infectados <- diff(infectados)
    
    d1_infectados <-
      ts(dados$Total,
         start = dados[1, 'Dia_Juliano'] - 1,
         frequency = 1)
    
    dif2_infectados <- diff(d1_infectados)
    
    tc <- dif2_infectados / dif_infectados - 1
    tc[is.nan(tc)] <- NA
    tc[is.infinite(tc)] <- NA
    
    # Preenchimento de falhas
    aux_kalman <- NULL
    aux_kalman <-
      try(tc <- na_kalman(tc, model = "StructTS", smooth = TRUE))
    
    # Posta resultados no Twitter.
    jpeg(
      "mm7_temp.jpg",
      width = 1400,
      height = 650,
      quality = 90
    )
    if (!is.null(aux_kalman)) {
      m7_tp <- sma(
        tc,
        order = 7,
        h = 1,
        interval = 'none',
        silent = FALSE
      )
      abline(h = 1,
             col = "gray",
             lty = 'dashed')
      dev.off()
      
      try(prev_m7_tp <- predict(m7_tp)$forecast[1])
      
      # Estatistica por hora não utilizada, muito volátil.
      
      #try(updateStatus(
      #  paste(
      #    format(aux_dh, '%d/%m/%Y'),
      #    '-',
      #    aux_nivel,
      #    '- Taxa de crescimento diario prevista (t+1) com SMA(7):',
      #    formatC(
      #      last(m7_tp$forecast),
      #      format = 'f',
      #      digits =
      #        5
      #    )
      #  ),
      #  mediaPath = "mm7_temp.jpg"
      #))
    }
    
    # Velocidade de avanco de infectados
    jpeg(
      "in1s_temp.jpg",
      width = 1400,
      height = 650,
      quality = 90
    )
    
    ts_in1s <- va(infectados)
    try(modelo_in1s <- sma(
      tail(ts_in1s, 14),
      order = 7,
      h = 1,
      interval = 'none',
      silent = FALSE
    ))
    # abline(h = 0, col = 'gray', lty = 'dashed')
    
    dev.off()
    
    try(updateStatus(
      paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Velocidade de avanço para casos novos (t+1) com SMA(7):',
        formatC(
          last(modelo_in1s$forecast),
          format = 'f',
          digits =  5
        )
      ),
      mediaPath = "in1s_temp.jpg"
    ))
    
    # Estima R
    dados_MT <- as.numeric(diff(infectados))
    
    try(res <- estimate_R(dados_MT, method = "parametric_si",
                          config = make_config(list(
                            mean_si = 3.96, std_si = 4.75
                          ))))
    
    jpeg(
      'r_temp.jpg',
      width = 700,
      height = 500,
      quality = 95
    )
    plot(res)
    dev.off()
    
    try({
      r_aux_025 <- tail(res$R$'Quantile.0.025(R)', 1)
      r_aux_975 <- tail(res$R$'Quantile.0.975(R)', 1)
      r_aux_50 <- tail(res$R$'Median(R)', 1)
      r_saida <- paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Valor de R[P02.5, P50, P97.5]: [',
        formatC(r_aux_025, format = 'f', digits = 2),
        ', ',
        formatC(r_aux_50,  format = 'f', digits = 2),
        ', ',
        formatC(r_aux_975, format = 'f', digits = 2),
        ']'
      )
    })
    
    try(updateStatus(r_saida, mediaPath = 'r_temp.jpg'))
    
    # Insere registro na tabela de Rs por ERS
    df_lolipop <-
      rbind(df_lolipop,
            data.frame(aux_nivel, r_aux_025, r_aux_975, r_aux_50))
    
    # Avalia a média móvel do R dos últimos 30 dias
    r30 <- tail(res$R$'Median(R)', 30)
    
    jpeg(
      "r30_temp.jpg",
      width = 1400,
      height = 650,
      quality = 90
    )
    
    mm_r30 <- sma(
      r30,
      order = 7,
      h = 1,
      interval = 'none',
      silent = "none"
    )
    dev.off()
    
    # Classificacao do R
    Rt <- last(mm_r30$forecast)
    if (Rt < 1.0) {
      risco_r <- 'Muito baixo'
    } else {
      if (Rt < 1.3) {
        risco_r <- 'Baixo'
      } else {
        if (Rt < 1.7) {
          risco_r <- 'Moderado'
        } else {
          if (Rt < 2.5) {
            risco_r <- 'Alto'
          } else {
            risco_r <- 'Muito alto'
          }
        }
      }
    }
    
    # Avalia a tendência com base nas medias moveis dos ultimos 7 dias.
    tendencia <-
      cor.test(time(tail(res$R$'Median(R)', 7)), tail(res$dates, 7), method = "spearman")
    res_tendencia <- ''
    if ((tendencia$p.value < 0.005) && (abs(tendencia$estimate) > 0.5)) {
      if (tendencia$estimate > 0) {
        res_tendencia <- '+'
      } else {
        res_tendencia <- '-'
      }
    }
    
    # Posta resultados no Twitter.
    try(updateStatus(paste(
      format(aux_dh, '%d/%m/%Y'),
      '-',
      aux_nivel,
      'Rt:',
      '- Risco:',
      risco_r,
      res_tendencia
    ),
    mediaPath = "r30_temp.jpg"))
    
    # Estima estatistica do Estado
    if (aux_nivel == 'Mato Grosso') {
      # Velocidade de avanco em enfermarias utilizadas
      enfermarias <-
        ts(
          total_param_mt$enfermaria_utilizadas,
          start = total_param_mt[1, 'Dia_Juliano'] - 1,
          frequency = 1
        )
      
      # Posta resultados no Twitter.
      jpeg(
        "enfermarias_temp.jpg",
        width = 1400,
        height = 650,
        quality = 90
      )
      
      ts_enfermarias <- va(enfermarias)
      modelo_enfermarias <- sma(
        tail(ts_enfermarias, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = FALSE
      )
      abline(h = 0, col = 'gray', lty = 'dashed')
      
      dev.off()
      
      try(updateStatus(
        paste(
          format(aux_dh, '%d/%m/%Y'),
          '-',
          aux_nivel,
          '- Velocidade de avanço para enfermarias em MT (t+1) com SMA(7):',
          formatC(
            last(modelo_enfermarias$forecast),
            format = 'f',
            digits =
              5
          )
        ),
        mediaPath = "enfermarias_temp.jpg"
      ))
      
      # Velocidade de avanco em UTIs
      UTI <-
        ts(total_param_mt$uti_utilizadas,
           start = total_param_mt[1, 'Dia_Juliano'] - 1,
           frequency = 1)
      
      jpeg(
        "UTI_temp.jpg",
        width = 1400,
        height = 650,
        quality = 90
      )
      
      ts_UTI <- va(UTI)
      modelo_UTI <- sma(
        tail(ts_UTI, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = FALSE
      )
      abline(h = 0, col = 'gray', lty = 'dashed')
      
      dev.off()
      
      try(updateStatus(
        paste(
          format(aux_dh, '%d/%m/%Y'),
          '-',
          aux_nivel,
          '- Velocidade de avanço para UTIs em MT (t+1) com SMA(7):',
          formatC(
            last(modelo_UTI$forecast),
            format = 'f',
            digits =
              5
          )
        ),
        mediaPath = "UTI_temp.jpg"
      ))
    }
  }
}

# Salva as projecoes por aglomerado
write.table(
  file = 'prev_aglomerados.csv',
  proj_aglomerado,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)


# Barras de erro
jpeg(
  'lolipop_temp.jpg',
  width = 1300,
  height = 800,
  quality = 95
)

# Grafico do R
try({
  ggplot(df_lolipop) +
    geom_segment(aes(
      x = aux_nivel,
      xend = aux_nivel,
      y = r_aux_025,
      yend = r_aux_975
    ),
    color = "dark blue") +
    geom_point(
      aes(x = aux_nivel, y = r_aux_025),
      color = rgb(0.2, 0.7, 0.1, 0.5),
      size = 3
    ) +
    geom_point(
      aes(x = aux_nivel, y = r_aux_975),
      color = rgb(0.7, 0.2, 0.1, 0.5),
      size = 3
    ) +
    geom_point(
      aes(x = aux_nivel, y = r_aux_50),
      color = 'black',
      size = 1,
      shape = 3
    ) +
    geom_hline(
      aes(yintercept = 1),
      color = 'black',
      size = 1,
      linetype = "dashed"
    ) +
    coord_flip() +
    labs(
      x = 'Escritórios Regionais de Saúde (ERS)',
      y = 'Número básico de reprodução (R)',
      title = "Número básico de reprodução por Escritorio Regional de Saúde em Mato Grosso",
      subtitle = "Estimativas a partir da série temporal",
      caption = paste('Processamento: ', format(aux_dh, '%d/%m/%Y'))
    ) +
    theme_ipsum() +
    theme(legend.position = "none", )
  dev.off()
})

# Posta intervalosde confianca do R no Twitter
r_saida <-
  paste(format(aux_dh, '%d/%m/%Y'),
        '- Intervalos de confiança do R por ERS e Mato Grosso.')
try(updateStatus(r_saida, mediaPath = 'lolipop_temp.jpg'))


#
# SIR
#

# Funcao SIR
sir <- function(suscetiveis,
                real_infectados,
                recuperados,
                dias_rec,
                aux_inf.prob,
                aux_act.rate,
                n_passos,
                n_simulacoes,
                n_cores = 3) {
  # Define os parametros de controle
  controle <- control.icm(
    type = "SIR",
    nsteps = n_passos,
    nsims = n_simulacoes,
    ncores = n_cores
  )
  
  # Inicializa a simulacao
  init <-
    init.icm(
      s.num = as.numeric(suscetiveis),
      i.num = as.numeric(real_infectados),
      r.num = as.numeric(recuperados)
    )
  
  # Define os parametros epidemiologicos
  param <-
    param.icm(
      inf.prob = aux_inf.prob,
      act.rate = aux_act.rate,
      
      rec.rate = (1 - tm / 100) / dias_rec,
      # Taxa media diaria de recuperacao com imunidade.
      
      a.rate = (1.5 * 6.08 / 1000) / 365,
      # Taxa de chegada ou entrada, 50% maior que a de entrada para considerar nascimentos e imigracao.
      
      ds.rate = (0.5 / 100) / 365,
      # Taxa de mortalidade diaria de suspeitos.
      
      di.rate = tm / 22,
      # Taxa de mortalidade diaria de infectados.
      
      dr.rate = (6.08 / 1000) / 365
    ) # Taxa de mortalidade diaria dos recuperados.
  
  sim <- icm(param, init, controle)
  return(sim)
}

# Configuracao da simulacao
n_passos <- 240
n_simulacoes <- 15
n_cores <- 16

df_sim <- data.frame(
  time = integer(),
  s.num = double(),
  i.num = double(),
  r.num = double(),
  tm = double(),
  dias_rec  = double(),
  disponiveis_atendimento = double(),
  inf.prob = double(),
  act.rate = double(),
  data_processamento = character(),
  data = as.Date(character()),
  aglomerado = character(),
  regiao = character(),
  taxa_ocupacao = double(),
  perc_atendimento = double(),
  perc_uti = double(),
  perc_mortalidade = double(),
  perc_nao_detectados = double(),
  infectados_oficial = double(),
  recuperados_oficial = double(),
  disponiveis_uti = double(),
  disponiveis_atendimento = double(),
  stringsAsFactors = FALSE
)

# Seta data de processamento e nome de arquivo com resultados
data_processamento <- format(Sys.time(), '%d%m%Y-%H%M%S')
aux_arquivo <- paste('aux_sim', data_processamento, '.csv')

for (aux_aglomerados in 1:length(aglomerados[['aglomerado']])) {
  # Seleciona registro do municipio/aglomerado
  aux_municipio <-
    aux_dados[aux_dados$Municipio == aglomerados[aux_aglomerados, 'aglomerado'],]
  aux_municipio <- aux_municipio[nrow(aux_municipio), ]
  
  # Seta infectados/recuperados+obitos/sucetiveis
  infectados <- aux_municipio$Infectados
  recuperados <- aux_municipio$Recuperados + aux_municipio$Obitos
  suscetiveis <-
    aglomerados[['suscetiveis']][aux_aglomerados] - infectados - recuperados
  # suscetiveis <- 300 # Para testar/acelerar processamento
  
  # Nome do aglomerado
  aglomerado <- aglomerados[['aglomerado']][aux_aglomerados]
  
  # Seta TM
  tm <- aglomerados[['tm']][aux_aglomerados] / 100
  
  # Define os parametros do SIR por aglomerado
  dias_rec <- aglomerados[['dias_rec']][aux_aglomerados]
  aux_inf.prob <- aglomerados[['inf_prob']][aux_aglomerados]
  aux_act.rate <- aglomerados[['act_rate']][aux_aglomerados]
  
  # Seta parâmetros do sistema de saude local
  leitos_uti <- aglomerados[['leitos_uti']][aux_aglomerados]
  leitos_atendimento <-
    aglomerados[['leitos_atendimento']][aux_aglomerados]
  taxa_ocupacao <- aglomerados[['taxa_ocupacao']][aux_aglomerados]
  perc_atendimento <-
    aglomerados[['perc_atendimento']][aux_aglomerados]
  perc_uti <- aglomerados[['perc_uti']][aux_aglomerados]
  perc_mortalidade <-
    aglomerados[['perc_mortalidade']][aux_aglomerados]
  pind_infectados <-
    aglomerados[['perc_nao_detectados']][aux_aglomerados]
  
  # Seta parametros de subnotificacao
  real_infectados <-
    round(infectados + infectados / (1 - pind_infectados), digits = 0)
  real_recuperados <-
    round(recuperados + recuperados / (1 - pind_infectados), digits = 0)
  
  # Processa a simulacao.
  sim <-
    sir(
      suscetiveis,
      real_infectados,
      real_recuperados,
      dias_rec,
      aux_inf.prob,
      aux_act.rate,
      n_passos,
      n_simulacoes,
      n_cores
    )
  
  # Indica que finalizou o processamento
  print(paste('SIR -', aglomerado))
  
  # Cria dataframe com resultados brutos das simulacoes.
  aux_df_sim <-
    as.data.frame(sim)[, c('time', 's.num', 'i.num', 'r.num')]
  aux_df_sim['inf.prob'] <- aux_inf.prob
  aux_df_sim['act.rate'] <- aux_act.rate
  aux_df_sim['dias_rec'] <- dias_rec
  aux_df_sim['tm'] <- tm
  
  # Agrega as simulacoes.
  aux_df_sim <-
    aggregate(aux_df_sim, list(aux_df_sim$time), mean)
  
  # Cria campo com data e hora de processamento.
  aux_df_sim['data_processamento'] <- data_processamento
  
  # Insere variaveis da simulacao.
  # Data da simulacao.
  aux_df_sim['data'] <-
    aux_df_sim$time + aglomerados[['data']][aux_aglomerados]
  # Nome do aglomerado/municipio.
  aux_df_sim['aglomerado'] <-
    aglomerados[['aglomerado']][aux_aglomerados]
  # Taxa de ocupacao.
  aux_df_sim['taxa_ocupacao'] <- taxa_ocupacao
  # Percentual de individuos que precisarao de atendimento hospitalar.
  aux_df_sim['perc_atendimento'] <- perc_atendimento
  # Percentual de individuos que precisarao de UTI.
  aux_df_sim['perc_uti'] <- perc_uti
  # Percentual de individuos que irao a obito.
  aux_df_sim['perc_mortalidade'] <- perc_mortalidade
  # Percentual de infectados assintomaticos.
  aux_df_sim['perc_nao_detectados'] <- pind_infectados
  # Numero de infectados oficial.
  aux_df_sim['infectados_oficial'] <-
    round((aux_df_sim['i.num'] - aux_df_sim['i.num'] * pind_infectados) / (2 - pind_infectados))
  # Numero de recuperados oficial.
  aux_df_sim['recuperados_oficial'] <-
    round((aux_df_sim['r.num'] - aux_df_sim['r.num'] * pind_infectados) / (2 - pind_infectados))
  # Numero de UTIs disponiveis.
  aux_df_sim['disponiveis_uti'] <-
    leitos_uti - aux_df_sim['i.num'] * perc_uti
  # Numero de enfermarias disponiveis.
  aux_df_sim['disponiveis_atendimento'] <-
    leitos_atendimento - aux_df_sim['i.num'] * perc_atendimento
  # ERS
  aux_df_sim['regiao'] <-
    aglomerados[['regiao']][aux_aglomerados]
  
  # Adiciona resultados em arquivo.
  aux_df_sim <- aux_df_sim[,-1]
  df_sim <- rbind(df_sim, aux_df_sim)
}

# Salva arquivo de arquivo de resultados
write.table(
  df_sim,
  aux_arquivo,
  append = TRUE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

# Agrega todos os casos por dia e o associa ao MT.
mt_df_sim <-
  aggregate(df_sim[, c('time',
                       's.num',
                       'infectados_oficial',
                       'recuperados_oficial')], list(df_sim$time), sum)[-2]
mt_df_sim$ERS <- 'Mato Grosso'
colnames(mt_df_sim) <-
  c('Data',
    'Suscetiveis',
    'Infectados',
    'Recuperados/Fatalidades',
    'ERS')
mt_df_sim$ERS <- factor(mt_df_sim$ERS)
mt_df_sim$Data <- aux_dh + mt_df_sim$Data - 1

# Agrega os casos por regiao (ERS) e dia.
aux_df_sim <-
  aggregate(df_sim[, c('time',
                       's.num',
                       'infectados_oficial',
                       'recuperados_oficial')],
            list(df_sim$regiao,
                 df_sim$time),
            sum)[-3]
colnames(aux_df_sim) <-
  c('ERS',
    'Data',
    'Suscetiveis',
    'Infectados',
    'Recuperados/Fatalidades')
aux_df_sim$ERS <- factor(aux_df_sim$ERS)
aux_df_sim$Data <- aux_dh + aux_df_sim$Data - 1

# Localiza ponto de maximo no Mato Grosso e o armazena
ers_aux_df_sim <- mt_df_sim
ers_aux_df_sim <-
  ers_aux_df_sim[which.max(ers_aux_df_sim$Infectados), ]

# Localiza os maximos de infectados nas ERS e os armazenam
for (aux_ers in levels(aux_df_sim$ERS)) {
  aux_ers_df_sim <- aux_df_sim[aux_df_sim$ERS == aux_ers,]
  aux_ers_df_sim <-
    aux_ers_df_sim[which.max(aux_ers_df_sim$Infectados), ]
  ers_aux_df_sim <- rbind(ers_aux_df_sim, aux_ers_df_sim)
}

# Cria post para o Twitter
jpeg(
  "pico_temp.jpg",
  width = 300,
  height = 750,
  quality = 90
)
grid.table(data.frame(ers_aux_df_sim[, c('Data', 'ERS')]))
dev.off()
r_saida <-
  paste(
    format(aux_dh, '%d/%m/%Y'),
    '- Estimativa das datas e valores de maximo de infectados por ERS e Mato Grosso.'
  )
try(updateStatus(r_saida, mediaPath = 'pico_temp.jpg'))

# Estima o total de UTIs utilizadas no MT
# Seleciona o ultimo registro do arquivo de parâmetros de MT.
param_mt <- tail(total_param_mt, 1)

# Define atalhos para os parâmetros.
uti_sus <- param_mt$leitos_uti_sus
p_uti <- param_mt$uti_utilizadas/param_mt$casos_confirmados
uti_utilizados <- param_mt$uti_utilizadas
enfermaria_utilizados <- param_mt$enfermaria_utilizadas
in_to_uti <- param_mt$to_leitos_uti_sus

# Cria post no Twitter.
jpeg(
  "cepi_temp.jpg",
  width = 1200,
  height = 650,
  quality = 90
)
aux_mt_df_sim <-
  mt_df_sim[1:(which.max(mt_df_sim$Infectados) + 30), ]
plot(
  aux_mt_df_sim$Data,
  aux_mt_df_sim$'Recuperados/Fatalidades' / 10000,
  type = "l",
  frame = TRUE,
  cex = 0.1,
  pch = 1,
  col = "blue",
  xlab = "Data",
  ylab = "Indivíduos/10'000",
  lty = 1,
  lwd = 1
)
lines(
  aux_mt_df_sim$Data,
  aux_mt_df_sim$Infectados / 10000,
  cex = 0.1,
  pch = 1,
  col = "red",
  type = "l",
  lty = 1,
  lwd = 1
)
legend(
  "topleft",
  legend = c("Recuperados/Fatalidades", "Infectados"),
  box.lty = 1,
  col = c("blue", "red"),
  lty = 1,
  cex = 0.8
)
dev.off()

r_saida <-
  paste(
    format(aux_dh, '%d/%m/%Y'),
    '- Projeção de infectados e recuperados/fatalidades no Mato Grosso.'
  )
try(updateStatus(r_saida, mediaPath = 'cepi_temp.jpg'))

# Cria indexador
mt_df_sim$d <-
  as.data.frame(as.numeric(mt_df_sim$Data - first(mt_df_sim$Data)))
colnames(mt_df_sim)[6] <- 'd'

# Funcao de decaimento exponencial adaptada para ajustar a ocupacao pre-existente.
decaimento <- function(d, p = 1.0508) {
  # Decaimento de 50% em 14 dias: 1.0508.
  # Decaimento de 99% em 20 dias: 1.257.
  decaimento <- 1 - (1 / p) ** d
  return(decaimento)
}
# Calcula a taxa de ocupacao das UTIs e cria um indice para o valor mais proximo.
# Por algum motivo o R esta entendendo que as colunas sao listas.
mt_df_sim$to_uti <- as.numeric(unlist(
  (
    uti_sus * in_to_uti * (1 - decaimento(mt_df_sim$d)) +
      p_uti * mt_df_sim$Infectados * decaimento(mt_df_sim$d)
  ) /
    uti_sus
))

data_max_mt_df_sim <-
  mt_df_sim[which.max(mt_df_sim$to_uti), 'Data']
mt_df_sim$dif_data <- sign(mt_df_sim$Data - data_max_mt_df_sim)

# Localiza o registro mais proximo da taxa de ocupacao de 60%.
mt_df_sim$abs_dif_to_uti <- as.numeric(unlist(mt_df_sim$dif_data *
                                                1 / abs(0.6 - mt_df_sim$to_uti)))
l60_mt_df_sim <- which.min (mt_df_sim$abs_dif_to_uti)
p60_mt_df_sim <- mt_df_sim[l60_mt_df_sim, ]

# Cria post no Twitter
jpeg('uti_temp.jpg',
     width = 1300,
     height = 800,
     quality = 95)
plot(
  mt_df_sim$Data,
  mt_df_sim$to_uti*100,
  type = "l",
  frame = TRUE,
  cex = 0.1,
  pch = 1,
  col = "blue",
  xlab = "Data",
  ylab = "TO de UTI/SUS/COVID-19 Mato Grosso (%)",
  lty = 1,
  lwd = 1
)
abline(h = 60, col = "red", lty = 'dashed')
dev.off()

r_saida <-
  paste(
    format(aux_dh, '%d/%m/%Y'),
    '- MT - Projeção da data em que a TO/SUS/COVID-19 atingirá 60% da capacidade atual',
    p60_mt_df_sim$Data
  )
try(updateStatus(r_saida, mediaPath = 'uti_temp.jpg'))

toc()


#
# Monte Carlo
#


tic()

n_simulacoes <-
  5 # Numero de simulacoes por iteração para o calculo de media.
n_cores <- 23 # Numero de cores disponiveis no computador.
mc_n_simulacoes <- 200 # Numero de simulacoes por MC.
max_dias <- 14 # Horizonte a ser considerado.
min_casos <- 5# Minimo de infectados para se efetuar a busca de parametros

for (aux_aglomerados in 1:length(aglomerados[['aglomerado']])) {
  aglomerado <- aglomerados[['aglomerado']][aux_aglomerados]
  serie_municipio <-
    aux_dados[aux_dados$Municipio == aglomerado, ]
  aux_municipio <-
    serie_municipio[nrow(aux_municipio),]
  
  # infectados <- aglomerados[['infectados']][aux_aglomerados]
  # recuperados <- aglomerados[['recuperados']][aux_aglomerados]
  
  # Seleciona serie
  #serie_municipio <-
  #  mt_dados[mt_dados$Municipio == as.character(aux_municipio$Municipio),]
  # serie_municipio <- total_aux_municipio
  
  # Seleciona o primeiro registro
  inicio <-
    serie_municipio[which.min(serie_municipio$Dia_Juliano),]
  
  # Seleciona o ultimo registro
  fim <- serie_municipio[which.max(serie_municipio$Dia_Juliano),]
  
  # Número de passos da simulacao
  n_passos <- fim$Dia_Juliano - inicio$Dia_Juliano + 1
  
  # Parametros do SIR por aglomerado
  dias_rec <- aglomerados[['dias_rec']][aux_aglomerados]
  aux_inf.prob <- aglomerados[['inf_prob']][aux_aglomerados]
  aux_act.rate <- aglomerados[['act_rate']][aux_aglomerados]
  
  # Critério de mínimo para iniciar a reavaliar os parametros
  if ((n_passos > max_dias) && (fim$Total > min_casos)) {
    print(paste('MC - ', aglomerado))
    
    # Parametros iniciais
    ref_serie_municipio <-
      serie_municipio[serie_municipio$Dia_Juliano > fim$Dia_Juliano - max_dias,]
    inicio <-
      ref_serie_municipio[which.min(ref_serie_municipio$Dia_Juliano),]
    n_passos <- fim$Dia_Juliano - inicio$Dia_Juliano + 1
    infectados <- inicio$Infectados
    recuperados <- inicio$Recuperados + inicio$Obitos
    suscetiveis <-
      aglomerados[['suscetiveis']][aux_aglomerados] - infectados - recuperados
    # suscetiveis <- 300 # Para teste
    
    tm <- aglomerados[['tm']][aux_aglomerados] / 100
    
    # Covariaveis
    leitos_uti <- aglomerados[['leitos_uti']][aux_aglomerados]
    leitos_atendimento <-
      aglomerados[['leitos_atendimento']][aux_aglomerados]
    taxa_ocupacao <-
      aglomerados[['taxa_ocupacao']][aux_aglomerados]
    perc_atendimento <-
      aglomerados[['perc_atendimento']][aux_aglomerados]
    perc_uti <- aglomerados[['perc_uti']][aux_aglomerados]
    perc_mortalidade <-
      aglomerados[['perc_mortalidade']][aux_aglomerados]
    pind_infectados <-
      aglomerados[['perc_nao_detectados']][aux_aglomerados]
    
    real_infectados <-
      infectados + infectados / (1 - pind_infectados)
    real_recuperados <-
      recuperados + recuperados / (1 - pind_infectados)
    
    min_mse_sir <- Inf
    for (aux_mc_n_simulacoes in 1:mc_n_simulacoes) {
      # Processa a simulacao
      sim <-
        sir(
          suscetiveis,
          real_infectados,
          real_recuperados,
          dias_rec,
          aux_inf.prob,
          aux_act.rate,
          n_passos,
          n_simulacoes,
          n_cores
        )
      
      print(paste('MC -', aux_mc_n_simulacoes, '-', aglomerado))
      aux_df_sim <-
        as.data.frame(sim)[, c('time', 's.num', 'i.num', 'r.num')]
      aux_df_sim['inf.prob'] <- aux_inf.prob
      aux_df_sim['act.rate'] <- aux_act.rate
      aux_df_sim['dias_rec'] <- dias_rec
      aux_df_sim['tm'] <- tm
      
      aux_df_sim <-
        aggregate(aux_df_sim, list(aux_df_sim$time), mean)
      
      aux_df_sim <- aux_df_sim[, -1]
      
      aux_df_sim['infectados_oficial'] <-
        round((aux_df_sim['i.num'] - aux_df_sim['i.num'] * pind_infectados) / (2 - pind_infectados))
      aux_df_sim['recuperados_oficial'] <-
        round((aux_df_sim['r.num'] - aux_df_sim['r.num'] * pind_infectados) / (2 - pind_infectados))
      
      aux14_df_sim <-
        aux_df_sim[aux_df_sim$time >= n_passos - max_dias, ]
      aux14_df_sim$time <-
        aux14_df_sim$time + inicio$Dia_Juliano - 1
      
      aux14_municipio <-
        aux_dados[(aux_dados$Dia_Juliano > fim$Dia_Juliano - max_dias) &
                    (aux_dados$Municipio == aglomerado), ]
      
      aux_14 <-
        merge(aux14_df_sim,
              aux14_municipio,
              by.x = 'time',
              by.y = 'Dia_Juliano')[, c(
                'infectados_oficial',
                'recuperados_oficial',
                'Infectados',
                'Recuperados',
                'Obitos'
              )]
      
      aux_14$mse <-
        (aux_14$infectados_oficial - aux_14$Infectados) ** 2 +
        (aux_14$recuperados_oficial - aux_14$Recuperados - aux_14$Obitos) ** 2
      
      mse <- mean(aux_14$mse)
      
      if (mse < min_mse_sir) {
        min_mse_sir <- mse
        atm_par_mse <-
          data.frame(aglomerado,
                     aux_inf.prob,
                     aux_act.rate,
                     dias_rec,
                     mse)
        colnames(atm_par_mse) <-
          c('municipio',
            'inf.prob',
            'act.rate',
            'dias_rec',
            'mse')
      }
      
      aux_inf.prob <-
        rnorm(n = 1,
              mean = aux_inf.prob,
              sd = aux_inf.prob / 15)
      dias_rec <- rnorm(n = 1,
                        mean = dias_rec,
                        sd = dias_rec / 15)
      aux_act.rate <-
        rnorm(n = 1,
              mean = aux_act.rate,
              sd = aux_act.rate / 15)
    }
    atm_par_mse$inf.prob ->
      aglomerados[aglomerados$aglomerado == atm_par_mse$municipio, 'inf_prob']
    atm_par_mse$act.rate ->
      aglomerados[aglomerados$aglomerado == atm_par_mse$municipio, 'act_rate']
    atm_par_mse$dias_rec ->
      aglomerados[aglomerados$aglomerado == atm_par_mse$municipio, 'dias_rec']
    atm_par_mse$mse ->
      aglomerados[aglomerados$aglomerado == atm_par_mse$municipio, 'mse']
    
    par_df_sim <- rbind(par_df_sim, atm_par_mse)
  }
}

write.table(
  file = 'par_municipios.csv',
  par_df_sim,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

write.table(
  file = 'aglomerados.csv',
  aglomerados,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

write.xlsx(aglomerados, 'mcr_aglomerados.xlsx')

tproc <- round(as.numeric(toc()), digits = 1)

r_saida <-
  paste(format(aux_dh, '%d/%m/%Y'),
        '- MMC ref. concluido com',
        round(tproc[2] - tproc[1], 1),
        's.')
try(updateStatus(r_saida))


#
# Mapa de risco por ERS (em elaboracao)
#

# install.packages('rgdal')
# install.packages("abjutils")

library(rgdal)
library(maptools)
if (!require(gpclib))
  install.packages("gpclib", type = "source")
gpclibPermit()
library(ggplot2)
library(abjutils)

shp <-
  readOGR(dsn = 'ERS_MUNICIPIOS_MT.shp', stringsAsFactors = F)
shp <-
  readOGR(dsn = 'ERS_SAUDE_MT_FINAL.shp', stringsAsFactors = F)

aux_shp <- broom::tidy(shp, region = 'Ers')
aux_shp$i_msa <- tolower(rm_accent(aux_shp$id))

aux_ers <- data.frame(shp)
aux_ers$ers_msa <- tolower(rm_accent(aux_ers$Ers))

df_lolipop$ers_msa <- tolower(rm_accent(df_lolipop$aux_nivel))

aglomerados$i_msa <- tolower(rm_accent(aglomerados$aglomerado))
aglomerados$ers_msa <- tolower(rm_accent(aglomerados$regiao))

aglomerados_lolipop <- merge(aglomerados, df_lolipop, all.x=TRUE, by=c('ers_msa'))
aglomerados_lolipop <- aglomerados_lolipop[, c('i_msa', 'ers_msa', 'r_aux_50')]
colnames(aglomerados_lolipop) <- c('i_msa', 'ers_msa', 'R')
aux_shp2 <- merge(aux_shp, aglomerados_lolipop, all.x=TRUE, by=c('i_msa'))

nomes_ers <-
  aggregate(cbind(Longitude, Latitude) ~ Ers,
            data = aux_shp,
            FUN = mean)

map <- ggplot() +
  geom_polygon(data = aux_shp2,
               aes(
                 x = long,
                 y = lat,
                 group = id,
                 fill = R
               ),
               colour = "black") +
  #geom_text(data = nomes_ers,
  #          aes(x = Longitude, y = Latitude,
  #              label = Ers),
  #          size = 3) +
  theme_void()
plot(map)

