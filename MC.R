# Automatizacao de processamentos relativos ao COVID-19
# Pesquisa de otimo por Monte Carlo
# Autor: Carlo Ralph De Musis

n_simulacoes <- 5 # Numero de simulacoes por iteracoes para o calculo de media.
n_cores <- 3 # Numero de cores disponiveis no computador.
mc_n_simulacoes <- 200 # Numero de simulacoes por MC.
max_dias <- 14 # Horizonte a ser considerado.
min_casos <- 5 # Minimo de infectados para se efetuar a busca de parametros

# Define diretorio de trabalho e carrega progresso de infectados
# setwd("D:/Mega/COVID-19")
setwd("~/simulacoes/COVID-19/")

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
library(stringi)

# Configuracao do Twitter
try(setup_twitter_oauth(
  consumer_key = "K8H2usTrUWKWrdUDhqZoMtYxT",
  access_token = "1254974121265500160-ZW0Gy3GLNe1xWsMEasMUzQ0G67rufK",
  consumer_secret = "egNsRyTB1qLJciGT6qweJjFhc2ZOXdU4t86Zfx9x0G9hD1rgUa",
  access_secret = "YKNvlxrp5VhgVhJ3zZpP5KFalZPkKzC4VtnKS56JUoxtb"
))

# Liga o cronometro
tic("Inicio...")

# Define e cria dataframe para grafico de lolipop
df_lolipop <- data.frame(
  ers = character(),
  p025 = double(),
  p50 = double(),
  p975 = double(),
  stringsAsFactors = FALSE
)

par_df_sim <- data.frame(
  municipio = character(),
  inf.prob = double(),
  act.rate = double(),
  dias_rec = double(),
  mse = double(),
  stringsAsFactors = FALSE
)

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

# Carrega base de dados da ERSs
aux_dados <- read.csv("evolucao_ers.CSV", stringsAsFactors = FALSE)

aux_dados$Data <- dmy(aux_dados$Data)
aux_dados$Dia_Juliano <- yday(aux_dados$Data)
aux_dados$Total <-
  aux_dados$Infectados + aux_dados$Recuperados + aux_dados$Obitos

aux_dados <- aux_dados[!is.na(aux_dados[, 'Infectados']), ]

aux_dados$ERS <- factor(aux_dados$ERS)
mt_dados <- aux_dados[, c('ERS', 'Dia_Juliano', 'Total')]
mt_aux_dados <- aggregate(Total ~ ERS + Dia_Juliano, mt_dados, sum)

mt_total_dados <- aggregate(Total ~ Dia_Juliano, mt_dados, sum)
mt_total_dados$ERS <- 'Mato Grosso'
mt_aux_dados <- rbind(mt_aux_dados, mt_total_dados)


#
# Funcao SIR
#

sir <- function(suscetiveis,
                real_infectados,
                recuperados,
                dias_rec,
                aux_inf.prob,
                aux_act.rate, 
                n_passos, 
                n_simulacoes,
                n_cores=3) {
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

# Carrega aglomerados
# aglomerados <- read.xlsx("aglomerados.xlsx", 1, header = TRUE)
aglomerados <-
  read.csv("aglomerados.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

aglomerados$data <-
  as.Date(parse_date_time(aglomerados[["data"]], '%d%m%y'), format = '%d%m%y')

colnames(aglomerados)[1] <- 'aglomerado'

for (aux_aglomerados in 1:length(aglomerados[['aglomerado']])) {
  aglomerado <- aglomerados[['aglomerado']][aux_aglomerados]
  serie_municipio <-
    aux_dados[aux_dados$Municipio == aglomerado, ]
  #aux_municipio <-
  #  serie_municipio[nrow(aux_municipio),]
  
  # Seleciona o primeiro registro
  inicio <-
    serie_municipio[which.min(serie_municipio$Dia_Juliano),]
  
  # Seleciona o ultimo registro
  fim <- serie_municipio[which.max(serie_municipio$Dia_Juliano),]
  
  # Numero de passos da simulacao
  n_passos <- fim$Dia_Juliano - inicio$Dia_Juliano + 1
  
  # Parametros do SIR por aglomerado
  dias_rec <- aglomerados[['dias_rec']][aux_aglomerados]
  aux_inf.prob <- aglomerados[['inf_prob']][aux_aglomerados]
  aux_act.rate <- aglomerados[['act_rate']][aux_aglomerados]
  
  # Criterio de minimo para iniciar a reavaliar os parametros
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
          c('municipio', 'inf.prob', 'act.rate', 'dias_rec', 'mse')
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

data_processamento <- format(Sys.time(), '%d%m%Y')
aux_arquivo <- paste('smc', data_processamento, 
                     stri_rand_strings(1, 23, '[a-zA-Z0-9]'),
                     '.csv')

# Salva arquivo de configuracao dde aglomerados
write.table(
  file=aux_arquivo,
  aglomerados,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ',')


write.table(
  file='par_municipios.csv',
  par_df_sim,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

tproc <- round(as.numeric(toc()), digits = 1)
r_saida <-
  paste('MMC concluido com',
        round(tproc[2] - tproc[1], 1),
        's.')
print(r_saida)
try(updateStatus(r_saida))
