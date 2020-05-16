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
# library(openxlsx)
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

# Define e cria dataframe para gráfico de lolipop
df_lolipop <- data.frame(
  ers = character(),
  p025 = double(),
  p50 = double(),
  p975 = double(),
  stringsAsFactors = FALSE
)

# Configuração do Twitter
#try(setup_twitter_oauth(
#  consumer_key = "K8H2usTrUWKWrdUDhqZoMtYxT",
#  access_token = "1254974121265500160-ZW0Gy3GLNe1xWsMEasMUzQ0G67rufK",
#  consumer_secret = "egNsRyTB1qLJciGT6qweJjFhc2ZOXdU4t86Zfx9x0G9hD1rgUa",
#  access_secret = "YKNvlxrp5VhgVhJ3zZpP5KFalZPkKzC4VtnKS56JUoxtb"
#))

# Define diretório de trabalho e carrega progressão de infectados
setwd("D:/OneDrive/Notebooks/Python/Corona virus/SARIMA")
aux_dados <- read.csv("evolucao_ers.csv", stringsAsFactors = FALSE)

aux_dados$Data <- dmy(aux_dados$Data)
aux_dados$Dia_Juliano <- yday(aux_dados$Data)
aux_dados$Total <-
  aux_dados$Infectados + aux_dados$Recuperados + aux_dados$Obitos

aux_dados <- aux_dados[!is.na(aux_dados[, 'Infectados']),]

aux_dados$ERS <- factor(aux_dados$ERS)
mt_dados <- aux_dados[, c('ERS', 'Dia_Juliano', 'Total')]
mt_aux_dados <- aggregate(Total ~ ERS + Dia_Juliano, mt_dados, sum)

mt_total_dados <- aggregate(Total ~ Dia_Juliano, mt_dados, sum)
mt_total_dados$ERS <- 'Mato Grosso'
mt_aux_dados <- rbind(mt_aux_dados, mt_total_dados)

# mt_aux_dados <- mt_aux_dados[, -2]
# colnames(mt_aux_dados) <- c('Dia_Juliano', 'Total')
# sum(aux_dados[,])

# Resgata data do último registro.
# aux_dh <- Sys.time()
aux_dh <- aux_dados[nrow(aux_dados), 'Data']

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
  
  if (max(infectados) >= 10 && length(infectados) >= 10) {
    print(paste('TS -', aux_nivel))
    # Janela dos dados (se desejar)
    ## infectados <- window(infectados, start=15, end=end(infectados)[1])
    
    # Pequisa transformação de Box-Cox
    l <- BoxCox.lambda(infectados)
    
    # Testa os modelos disponíveis no pacote Forecast
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
      # nnetar_ajuste <- nnetar(infectados, repeats = 10, lambda = l)
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
      # nnetar_ajuste, # O modelo por RN apresenta alguns porblemas para validacao
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
    
    # Previsão para 3 dias
    previsao <- forecast(modelo, level = c(95), h = 3)
    
    jpeg(
      "grafico_temp.jpg",
      width = 1400,
      height = 700,
      quality = 90
    )
    
    #autoplot(
    #  previsao,
    #  ts.colour = 'dodgerblue4',
    #  predict.colour = 'blue',
    #  predict.linetype = 'dashed'
    #) +
    #  ggtitle(paste('Total de infectados para o Mato Grosso (', modelo$method, ')')) +
    #  xlab("Dias após o primeiro infectado") + ylab("Número de indivíduos") +
    #  theme_bw() +
    #  theme(plot.title = element_text(hjust = 0.5)) +
    #  scale_x_continuous(breaks = seq(
    #    from = 1,
    #    to =  end(infectados)[1] + 3,
    #    by = 1
    #  )) +
    #  scale_y_continuous(breaks = seq(
    #    from = 0,
    #    to =  previsao$upper[3],
    #    by = 10
    #  )) +
    #  scale_y_log10()
    
    plot(previsao,
         main = aux_nivel,
         xlab = "Dias",
         ylab = "Infectados")
    dev.off()
    
    
    # while (!is.null(dev.list())) dev.off()
    
    # Posta tweets da previsão de curto prazo
    
    try(updateStatus(
      paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Previsões de curto prazo  [Gráfico - 3 dias].'
      ),
      mediaPath = "grafico_temp.jpg"
    ))
    
    try(previsao <- forecast(modelo, level = c(95), h = 15))
    jpeg(
      "tabela_temp.jpg",
      width = 400,
      height = 450,
      quality = 90
    )
    grid.table(data.frame(previsao))
    dev.off()
    
    # Posta a tabela
    try({
      updateStatus(
        paste(
          format(aux_dh, '%d/%m/%Y'),
          '-',
          aux_nivel,
          '- Previsões de curto prazo [Tabela - 15 dias]'
        ),
        mediaPath = "tabela_temp.jpg"
      )
    })
    
    # Taxa de propagação
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
        h = 0,
        interval = 'none',
        silent = FALSE
      )
      dev.off()
      
      try(prev_m7_tp <- predict(m7_tp)$forecast[1])
      
      try(updateStatus(
        paste(
          format(aux_dh, '%d/%m/%Y'),
          '-',
          aux_nivel,
          '- Taxa de crescimento prevista (t+1), SMA(7):',
          formatC(
            last(m7_tp$fitted),
            format = 'f',
            digits =
              5
          )
        ),
        mediaPath = "mm7_temp.jpg"
      ))
    }
    
    # Estima R e posta tweets
    dados_MT <- as.numeric(diff(infectados))
    
    try({
      res <- estimate_R(dados_MT, method = "parametric_si",
                        config = make_config(list(
                          mean_si = 1.1, std_si = sd(dados_MT)
                        )))
      
      jpeg(
        'r_temp.jpg',
        width = 700,
        height = 500,
        quality = 95
      )
      plot(res)
      dev.off()
    })
    
    try({
      r_aux_025 <- tail(res$R$'Quantile.0.025(R)', 1)
      r_aux_975 <- tail(res$R$'Quantile.0.975(R)', 1)
      r_aux_50 <- tail(res$R$'Median(R)', 1)
      r_saida <- paste(
        format(aux_dh, '%d/%m/%Y'),
        '-',
        aux_nivel,
        '- Valor de R0 [P02.5, P50, P97.5]: [',
        formatC(r_aux_025, format = 'f', digits = 2),
        ', ',
        formatC(r_aux_50,  format = 'f', digits = 2),
        ', ',
        formatC(r_aux_975, format = 'f', digits = 2),
        ']'
      )
    })
    
    try(updateStatus(r_saida, mediaPath = 'r_temp.jpg'))
    
    df_lolipop <-
      rbind(df_lolipop,
            data.frame(aux_nivel, r_aux_025, r_aux_975, r_aux_50))
  }
}

# Seleciona a serie dos ultimos 30 dias da mediana do R (apenas nas sextas-feiras)
if (weekdays(aux_dh) == 'sexta-feira' && aux_nivel=="Mato Grosso") {
  r_tc_50 <- tail(res$R$'Median(R)', 30)
  
  jpeg(
    "mm7_temp.jpg",
    width = 1400,
    height = 650,
    quality = 90
  )
  
  m7_tp <- sma(
    r_tc_50,
    order = 7,
    h = 0,
    interval = 'none',
    silent = "none"
  )
  dev.off()

  # Classificação do R  
  if (last(m7_tp$fitted) < 1.0) {
    risco_r <- 'Muito baixo'
  } else {
    if (last(m7_tp$fitted) < 1.3) {
      risco_r <- 'Baixo'
    } else {
      if (last(m7_tp$fitted) < 1.7) {
        risco_r <- 'Moderado'
      } else {
        if (last(m7_tp$fitted) < 2.5) {
          risco_r <- 'Alto'
        } else {
          risco_r <- 'Muito alto'
        }
      }
    }
  }
  
  tendencia <- cor.test(time(tail(m7_tp$fitted, 7)), tail(m7_tp$fitted, 7), method="spearman")
  
  if (tendencia$p.value > 0.01) {
    res_tendencia <- ' '
  } else {
    if (tendencia$p.value > 0.005) {
      if (tendencia$estimate > 0) {
        res_tendencia <- '+'
      } else {
        res_tendencia <- '-'
      }
    } else {
      if (tendencia$estimate > 0) {
        res_tendencia <- '++'
      } else {
        res_tendencia <- '--'
      }
    }
}
  
  # Posta no Twitter
  try(updateStatus(
    paste(
      format(aux_dh, '%d/%m/%Y'),
      '-',
      aux_nivel,
      '- Série de R nos últimos 30 dias, SMA(7):',
      formatC(last(m7_tp$fitted),
              format = 'f',
              digits =
                5),
      ', risco:', risco_r, res_tendencia 
    ),
    mediaPath = "mm7_temp.jpg"
  ))
}

# Barras de erro
jpeg(
  'lolipop_temp.jpg',
  width = 1300,
  height = 800,
  quality = 95
)

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
      y = 'Número básico de reprodução (R0)',
      title = "Número básico de reprodução por Escritório Regional de Saúde em Mato Grosso",
      subtitle = "Estimativas a partir da série temporal",
      caption = paste('Processamento: ', format(aux_dh, '%d/%m/%Y'))
    ) +
    theme_ipsum() +
    theme(legend.position = "none",)
})

dev.off()

r_saida <-
  paste(format(aux_dh, '%d/%m/%Y'),
        '- Intervalos de confiança do R0 por ERS')
try(updateStatus(r_saida, mediaPath = 'lolipop_temp.jpg'))


# SIR

# Carrega aglomerados
# aglomerados <- read.xlsx("aglomerados.xlsx", 1, header = TRUE)
aglomerados <-
  read.csv("aglomerados.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

aglomerados$data <-
  as.Date(parse_date_time(aglomerados[["data"]], '%d%m%y'), format = '%d%m%y')

# Para resolver problemas de conversao
# aglomerados[2:6] <- lapply(aglomerados[2:6], as.numeric)

colnames(aglomerados)[1] <- 'aglomerado'

# COnfiguração da simulacao
n_passos <- 180
n_simulacoes <- 11
n_cores <- 3

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

data_processamento <- format(Sys.time(), '%d%m%Y-%H%M%S')
aux_arquivo <- paste('sim_mc', data_processamento, '.csv')

for (aux_aglomerados in 1:length(aglomerados[['aglomerado']])) {
  aux_municipio <-
    aux_dados[aux_dados$Municipio == aglomerados[aux_aglomerados, 'aglomerado'], ]
  aux_municipio <- aux_municipio[nrow(aux_municipio),]
  
  # infectados <- aglomerados[['infectados']][aux_aglomerados]
  # recuperados <- aglomerados[['recuperados']][aux_aglomerados]
  
  infectados <- aux_municipio$Infectados
  recuperados <- aux_municipio$Recuperados + aux_municipio$Obitos
  suscetiveis <-
    aglomerados[['suscetiveis']][aux_aglomerados] - infectados - recuperados
  # suscetiveis <- 300 # Para teste
  
  aglomerado <- aglomerados[['aglomerado']][aux_aglomerados]
  tm <- aglomerados[['tm']][aux_aglomerados] / 100
  
  # Parâmetros do sir por aglomerado
  dias_rec <- aglomerados[['dias_rec']][aux_aglomerados]
  aux_inf.prob <- aglomerados[['inf_prob']][aux_aglomerados]
  aux_act.rate <- aglomerados[['act_rate']][aux_aglomerados]
  
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
  
  real_infectados <-
    infectados + infectados / (1 - pind_infectados)
  
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
  
  aux_sim <- NULL
  aux_sim <- try(sim <- icm(param, init, controle))
  
  if (!is.null(aux_sim)) {
    print(paste('SIR -', aglomerado))
    aux_df_sim <-
      as.data.frame(sim)[, c('time', 's.num', 'i.num', 'r.num')]
    aux_df_sim['inf.prob'] <- aux_inf.prob
    aux_df_sim['act.rate'] <- aux_act.rate
    aux_df_sim['dias_rec'] <- dias_rec
    aux_df_sim['tm'] <- tm
    
    aux_df_sim <-
      aggregate(aux_df_sim, list(aux_df_sim$time), mean)
    
    aux_df_sim['data_processamento'] <- data_processamento
    aux_df_sim['data'] <-
      aux_df_sim$time + aglomerados[['data']][aux_aglomerados]
    aux_df_sim['aglomerado'] <-
      aglomerados[['aglomerado']][aux_aglomerados]
    aux_df_sim['taxa_ocupacao'] <- taxa_ocupacao
    aux_df_sim['perc_atendimento'] <- perc_atendimento
    aux_df_sim['perc_uti'] <- perc_uti
    aux_df_sim['perc_mortalidade'] <- perc_mortalidade
    aux_df_sim['perc_nao_detectados'] <- pind_infectados
    aux_df_sim['infectados_oficial'] <-
      (aux_df_sim['i.num'] - aux_df_sim['i.num'] * pind_infectados) / (2 - pind_infectados)
    aux_df_sim['recuperados_oficial'] <-
      (aux_df_sim['r.num'] - aux_df_sim['r.num'] * pind_infectados) / (2 - pind_infectados)
    aux_df_sim['disponiveis_uti'] <-
      leitos_uti - aux_df_sim['i.num'] * perc_uti
    aux_df_sim['disponiveis_atendimento'] <-
      leitos_atendimento - aux_df_sim['i.num'] * perc_atendimento
    
    aux_df_sim['regiao'] <-
      aglomerados[['regiao']][aux_aglomerados]
    
    aux_df_sim <- aux_df_sim[, -1]
    df_sim <- rbind(df_sim, aux_df_sim)
  }
}

write.table(
  df_sim,
  aux_arquivo,
  append = TRUE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

mt_df_sim <-
  aggregate(df_sim[, c('time', 's.num', 'i.num', 'r.num')], list(df_sim$time), sum)[-2]
mt_df_sim$ERS <- 'Mato Grosso'
colnames(mt_df_sim) <-
  c('Data',
    'Suscetiveis',
    'Infectados',
    'Recuperados/Fatalidades',
    'ERS')

aux_df_sim <-
  aggregate(df_sim[, c('time', 's.num', 'i.num', 'r.num')],
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
mt_df_sim$ERS <- factor(mt_df_sim$ERS)
aux_df_sim$Data <- aux_dh + aux_df_sim$Data - 1
mt_df_sim$Data <- aux_dh + mt_df_sim$Data - 1

ers_aux_df_sim <- mt_df_sim
ers_aux_df_sim <-
  ers_aux_df_sim[which.max(ers_aux_df_sim$Infectados), ]

for (aux_ers in levels(aux_df_sim$ERS)) {
  aux_ers_df_sim <- aux_df_sim[aux_df_sim$ERS == aux_ers,]
  aux_ers_df_sim <-
    aux_ers_df_sim[which.max(aux_ers_df_sim$Infectados), ]
  ers_aux_df_sim <- rbind(ers_aux_df_sim, aux_ers_df_sim)
}

jpeg(
  "pico_temp.jpg",
  width = 600,
  height = 350,
  quality = 90
)
grid.table(data.frame(ers_aux_df_sim))
dev.off()

r_saida <-
  paste(
    format(aux_dh, '%d/%m/%Y'),
    '- Estimativa das datas e valores de máximo de infectados por ERS'
  )
try(updateStatus(r_saida, mediaPath = 'pico_temp.jpg'))

toc()