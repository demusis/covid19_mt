# Analice de risco da COVID-19
# Carlo Ralph De Musis
# Junho/2020

import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['figure.figsize'] = [10, 7]

# Variavel antecedente: classificacao epidemiologica.
ce = ctrl.Antecedent(np.arange(0, 1, 0.05), 'Classificacao epidemiologica')

ce['Local'] = fuzz.trimf(ce.universe, [0, 0, 0.5])
ce['Indeterminada'] = fuzz.trimf(ce.universe, [0.25, 0.5, 0.75])
ce['Comunitaria']  = fuzz.trimf(ce.universe, [0.5, 1, 1])

# ce.view()

# Variavel antecedente: velocidade de avanco de casos novos.
vacn = ctrl.Antecedent(np.arange(-1, 2, 0.05), 'Velocidade de avanco de casos novos')

vacn['Desaceleracao'] = fuzz.trapmf(vacn.universe, [-1, -1, -0.5, 0])
vacn['Indeterminada'] = fuzz.trimf(vacn.universe, [-0.5, 0, 0.25])
vacn['Media']  = fuzz.trimf(vacn.universe, [0, 0.5, 1])
vacn['Alta']  = fuzz.trimf(vacn.universe, [0.75, 1, 1.5])
vacn['Muito alta']  = fuzz.trapmf(vacn.universe, [1, 1.5, 2, 2])

# vacn.view()


# Variavel antecedente: media movel do numero de reproducao basico.
rt = ctrl.Antecedent(np.arange(0, 3, 0.05), 'Media movel do numero basico de reproducao')

rt['Desaceleracao'] = fuzz.trapmf(rt.universe, [0, 0, 0.5, 1.0])
rt['Indeterminada'] = fuzz.trimf(rt.universe, [0.5, 1, 1.5])
rt['Media']  = fuzz.trimf(rt.universe, [1.25, 1.5, 1.75])
rt['Alta']  = fuzz.trimf(rt.universe, [1.5, 2.0, 2.5])
rt['Muito alta']  = fuzz.trapmf(rt.universe, [2.0, 2.5, 3, 3])

# rt.view()

# Variavel antecedente: possui infectados ativos.
ia = ctrl.Antecedent(np.arange(0, 1, 0.05), 'Possui infectados ativos')

ia['Nao'] = fuzz.trimf(ia.universe, [0, 0, 0.5])
ia['Indeterminado'] = fuzz.trimf(ia.universe, [0.25, 0.5, 0.75])
ia['Sim']  = fuzz.trimf(ia.universe, [0.5, 1, 1])

# ia.view()

# Variavel antecedente: taxa de ocupacao de UTI.
to_uti = ctrl.Antecedent(np.arange(0, 1, 0.05), 'Taxa de ocupacao de UTI')

to_uti['Muito baixa'] = fuzz.trapmf(to_uti.universe, [0, 0, 0.2, 0.4])
to_uti['Baixa'] = fuzz.trimf(to_uti.universe, [0.2, 0.4, 0.6])
to_uti['Media']  = fuzz.trimf(to_uti.universe, [0.4, 0.6, 0.75])
to_uti['Alta']  = fuzz.trimf(to_uti.universe, [0.6, 0.75, 0.9])
to_uti['Muito alta']  = fuzz.trapmf(to_uti.universe, [0.75, 0.9, 1, 1])

# to_uti.view()


# Risco
risco = ctrl.Consequent(np.arange(-17, 117, 1), 'Risco')

risco['Muito baixo'] = fuzz.trapmf(risco.universe, [-17, -17, 10, 20])
risco['Baixo'] = fuzz.gaussmf(risco.universe, 25, 5)
risco['Medio'] = fuzz.gaussmf(risco.universe, 50, 5)
risco['Alto'] = fuzz.gaussmf(risco.universe, 75, 5)
risco['Muito alto']  = fuzz.trapmf(risco.universe, [80, 90, 117, 117])

# risco.view()


#
# Regras
#

# Classificacao epidemiologica
regra_ce_1 = ctrl.Rule(ce['Local'] & ~ia['Nao'] 
                                   & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Baixo'])

# Infectados ativos.
regra_ia_1 = ctrl.Rule(ia['Nao'] & ~ce['Local'] 
                                 & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Baixo'])

# Transmissao e infectados ativos
regra_ce_ia_1 = ctrl.Rule(ce['Indeterminada'] & ia['Indeterminado'], risco['Medio'])
regra_ce_ia_2 = ctrl.Rule(ce['Local'] & ia['Nao'] 
                                      & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Muito baixo'])

# Velocidade de avanco
regra_vacn_1a = ctrl.Rule(vacn['Desaceleracao'] & ~(to_uti['Alta'] | to_uti['Muito alta'])
                                                & ~ce['Comunitaria'], risco['Muito baixo'])
regra_vacn_1b = ctrl.Rule(vacn['Desaceleracao'] & ~(to_uti['Alta'] | to_uti['Muito alta'])
                                                &  ce['Comunitaria'], risco['Alto'])
regra_vacn_1c = ctrl.Rule(vacn['Desaceleracao'] & (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])

regra_vacn_2a = ctrl.Rule(vacn['Indeterminada'] & ~(to_uti['Alta'] | to_uti['Muito alta'])
                                                & ~ce['Comunitaria'], risco['Baixo'])
regra_vacn_2b = ctrl.Rule(vacn['Indeterminada'] &  (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
regra_vacn_3a = ctrl.Rule(vacn['Media'] & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Medio'])
regra_vacn_3b = ctrl.Rule(vacn['Media'] &  (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
regra_vacn_4a = ctrl.Rule(vacn['Alta'] & ~to_uti['Muito alta'], risco['Alto'])
regra_vacn_4b = ctrl.Rule(vacn['Alta'] &  to_uti['Muito alta'], risco['Muito alto'])
regra_vacn_5 = ctrl.Rule(vacn['Muito alta'], risco['Muito alto'])

# Rt
regra_rt_1a = ctrl.Rule(rt['Desaceleracao'] & ~(to_uti['Alta'] | to_uti['Muito alta']) 
                                            & ~ce['Comunitaria'], risco['Muito baixo'])
regra_rt_1b = ctrl.Rule(rt['Desaceleracao'] & ~(to_uti['Alta'] | to_uti['Muito alta']) 
                                            &  ce['Comunitaria'], risco['Medio'])
regra_rt_1c = ctrl.Rule(rt['Desaceleracao'] & (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])

regra_rt_2a = ctrl.Rule(rt['Indeterminada'] & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Baixo'])
regra_rt_2b = ctrl.Rule(rt['Indeterminada'] &  (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
regra_rt_3a = ctrl.Rule(rt['Media'] & ~(to_uti['Alta'] | to_uti['Muito alta']), risco['Medio'])
regra_rt_3b = ctrl.Rule(rt['Media'] &  (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
regra_rt_4a = ctrl.Rule(rt['Alta'] & ~to_uti['Muito alta'], risco['Alto'])
regra_rt_4b = ctrl.Rule(rt['Alta'] &  to_uti['Muito alta'], risco['Muito alto'])
regra_rt_5 = ctrl.Rule(rt['Muito alta'], risco['Muito alto'])

# Velocidade de avanco e Rt
regra_vacn_rt_1 = ctrl.Rule(vacn['Muito alta']  & rt['Muito alta'] , risco['Muito alto'])
regra_vacn_rt_2 = ctrl.Rule(vacn['Alta']  & rt['Alta'] , risco['Muito alto'])
regra_vacn_rt_3 = ctrl.Rule(vacn['Media']  & rt['Media'] , risco['Alto'])
regra_vacn_rt_4a = ctrl.Rule(vacn['Indeterminada']  & rt['Indeterminada'] 
                                                    & ~(to_uti['Alta'] | to_uti['Muito alta']) 
                                                    & ~ce['Comunitaria'], risco['Medio'])
regra_vacn_rt_4b = ctrl.Rule(vacn['Indeterminada']  & rt['Indeterminada'] 
                                                    & (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
regra_vacn_rt_5a = ctrl.Rule(vacn['Desaceleracao']  & rt['Desaceleracao'] 
                                                    & ~(to_uti['Alta'] | to_uti['Muito alta']) 
                                                    & ~ce['Comunitaria'], risco['Muito baixo'])
regra_vacn_rt_5b = ctrl.Rule(vacn['Desaceleracao']  & rt['Desaceleracao'] 
                                                    & (to_uti['Alta'] | to_uti['Muito alta']), risco['Alto'])
# Taxa de ocupacao de UTI.
regra_to_uti_1 = ctrl.Rule(to_uti['Muito baixa'], risco['Muito baixo'])
regra_to_uti_2 = ctrl.Rule(to_uti['Baixa'], risco['Baixo'])
regra_to_uti_3 = ctrl.Rule(to_uti['Media'], risco['Medio'])
regra_to_uti_4 = ctrl.Rule(to_uti['Alta'], risco['Alto'])
regra_to_uti_5 = ctrl.Rule(to_uti['Muito alta'], risco['Muito alto'])



# Cria um sistema de controle e uma simulacao
risco_covid_ctrl = ctrl.ControlSystem([regra_ce_1, #ce
                                       
                                       regra_ia_1, # ia
                                       
                                       regra_ce_ia_1, regra_ce_ia_2, # ia e ce
                                       
                                       regra_vacn_1a, regra_vacn_1b, regra_vacn_1c, regra_vacn_2a, # vacn e to_uti
                                       regra_vacn_2b, regra_vacn_3a, regra_vacn_3b, regra_vacn_4a, 
                                       regra_vacn_4b, regra_vacn_5,
                                       
                                       regra_rt_1a, regra_rt_1b, regra_rt_1c, regra_rt_2a, #rt e to_uti
                                       regra_rt_2b, regra_rt_3a, regra_rt_3b, regra_rt_4a, 
                                       regra_rt_4b, regra_rt_5,
                                       
                                       regra_vacn_rt_1, regra_vacn_rt_2, regra_vacn_rt_3, regra_vacn_rt_4a, # vacn e rt
                                       regra_vacn_rt_4b, regra_vacn_rt_5a, regra_vacn_rt_5b,
                                       
                                       regra_to_uti_1, regra_to_uti_2, regra_to_uti_3, # to_uti 
                                       regra_to_uti_4, regra_to_uti_5])

risco_covid_sim = ctrl.ControlSystemSimulation(risco_covid_ctrl)



# Grafo das regras.
# risco_covid_ctrl.view()


# Grafo das variaveis linguísticas.
# risco_covid_ctrl.view_n()


def calculaRisco(aux):
    # Entrada.
    risco_covid_sim.input['Classificacao epidemiologica'] = aux[1]
    risco_covid_sim.input['Possui infectados ativos'] = aux[2]
    risco_covid_sim.input['Media movel do numero basico de reproducao'] = aux[3]
    risco_covid_sim.input['Velocidade de avanco de casos novos'] = aux[4]
    risco_covid_sim.input['Taxa de ocupacao de UTI'] = aux[5] 

    # Calculo.
    risco_covid_sim.compute()
    # print(risco_covid_sim.input)
    aux_risco = risco_covid_sim.output['Risco']
    if aux_risco < 0:
        aux_risco = 0
    elif aux_risco >100:
        aux_risco = 100    
    return(aux_risco)



# calculaRisco(['Teste', 0, -1, 0, 0, 0])



# Variaveis antecedentes.
# vacn.view(sim=risco_covid_sim)
# rt.view(sim=risco_covid_sim)
# to_uti.view(sim=risco_covid_sim)
# ia.view(sim=risco_covid_sim)

# Variavel consequente.
# risco.view(sim=risco_covid_sim)



df = pd.read_csv('.//resultados//2020-06-15 - bi_aglomerados.csv', decimal=',', sep=';')
# df.head()



# df = df[['municipio', 'classificacao_epidemiologica', 'infectados_n', 'r_mm7', 'p_va', 'to_uti']]
# df.head()



df.loc[~((df.classificacao_epidemiologica=='TRANSMISSAO COMUNITARIA') |
         (df.classificacao_epidemiologica=='TRANSMISSAO LOCAL'))
       , 'classificacao_epidemiologica'] = 0.5
df.loc[df.classificacao_epidemiologica=='TRANSMISSAO LOCAL', 'classificacao_epidemiologica'] = 0
df.loc[df.classificacao_epidemiologica=='TRANSMISSAO COMUNITARIA', 'classificacao_epidemiologica'] = 1


df.loc[np.isnan(df.infectados_n), 'infectados_n'] = 0.5
df.loc[df.infectados_n > 0, 'infectados_n'] = 1
df.loc[df.infectados_n == 0, 'infectados_n'] = 0

df.loc[np.isnan(df.r_mm7), 'r_mm7'] = 1
df.loc[np.isnan(df.p_va), 'p_va'] = 0

# df.head()


df['risco'] = df.apply(calculaRisco, axis = 1)

df.to_csv('risco_fuzzy.csv', decimal=',', sep=';')
df.to_excel('risco_fuzzy.xlsx')

