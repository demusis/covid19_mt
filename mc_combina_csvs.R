# Agrega arquivos .csv e os combina em arquivo de medias.

# install.packages('twitteR')
# install.packages('openxlsx')
 
library(openxlsx)
library(twitteR)

# Configuracao do Twitter
try(setup_twitter_oauth(
   consumer_key = "K8H2usTrUWKWrdUDhqZoMtYxT",
   access_token = "1254974121265500160-ZW0Gy3GLNe1xWsMEasMUzQ0G67rufK",
   consumer_secret = "egNsRyTB1qLJciGT6qweJjFhc2ZOXdU4t86Zfx9x0G9hD1rgUa",
   access_secret = "YKNvlxrp5VhgVhJ3zZpP5KFalZPkKzC4VtnKS56JUoxtb"
))

# Define diretorio de trabalho e carrega progresso de infectados
# setwd("D:/Mega/COVID-19")
# setwd("D:/Mega/COVID-19/Auxiliar")
setwd("~/simulacoes/COVID-19/")

nomes_csvs <- list.files(path = '.', pattern = glob2rx('mmc*.csv'))
datalist = lapply(nomes_csvs, function(x) {
  read.csv(
    file = x,
    header = TRUE,
    sep = ',',
    dec = '.',
    stringsAsFactors = FALSE
  )
})
df_sim <- Reduce(function(x, y) {
  rbind(x, y)
}, datalist)
data_processamento <- df_sim[1, 'data_processamento']

# Salva arquivo com todos os registros (bak)
write.table(
  df_sim,
  'dados_sim_mc.csv',
  append = FALSE,
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ','
)

df_sim$aglomerado <- factor(df_sim$aglomerado)

# Remove os arquivos temporarios
file.remove(list.files(path='.', pattern = glob2rx('sim_mc*.csv')))

aux_df_sim <- data.frame(df_sim[0, ])

for (aux_municipio in levels(df_sim$aglomerado)) {
  mun_df_sim <- df_sim[df_sim$aglomerado == aux_municipio,]
  aux_df_sim <- rbind(aux_df_sim, mun_df_sim[which.min(mun_df_sim$mse),])
}

write.xlsx(aux_df_sim, 'mc_aglomerados.xlsx')

write.table(
  file='mc_aglomerados.CSV',
  aux_df_sim,
  append = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = ',')

