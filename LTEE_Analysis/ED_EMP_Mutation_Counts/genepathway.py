import pandas as pd

# Import CSV of list of mutations from one clone/strain
df = pd.read_csv ('input_ara-2bMOTs.csv')

# Create output df
d = {'Generation': [500, 1000, 1500, 2000, 5000, 10000, 15000, 20000, 30000, 40000, 50000],
    'EMP': [0 for i in range(11)],
    'ED': [0 for i in range(11)],
    'OxPPP': [0 for i in range(11)],
    'NonOxPPP': [0 for i in range(11)],
    'Lower Glycolysis': [0 for i in range(11)],
    'Lower Glycolysis no pykF' : [0 for i in range(11)],
    'Glucose Transport': [0 for i in range(11)],
    'Met Biosynth': [0 for i in range(11)],
    'Homoserine Biosynth': [0 for i in range(11)],
    'Met + Homoser Biosynth': [0 for i in range(11)]}
df_output = pd.DataFrame(data=d)

# Create sets of genes
emp_set = {'pgi', 'pfkA', 'pfkB', 'fbaA', 'fbaB', 'tpiA'}
ed_set = {'edd', 'eda', 'pgl', 'zwf'}
oxppp_set = {'zwf', 'pgl', 'gndA'}
noxppp_set = {'rpe', 'rpiA', 'rpiB', 'tktA', 'tktB', 'talA', 'talB'}
lowgly_set = {'gapA', 'pgk', 'gpmA', 'gpmM', 'eno', 'pykA', 'ppsA', 'pykF'} 
lgnopyfk_set = {'gapA', 'pgk', 'gpmA', 'gpmM', 'eno', 'pykA', 'ppsA'}
glctrans_set = {'ptsG', 'ptsI', 'ptsP'}
met_set = {'metA', 'metB', 'metC', 'malY', 'metE', 'metH'}
hser_set = {'lysC', 'thrA', 'metL', 'asd'}
met_hser_set = {'metA', 'metB', 'metC', 'malY', 'metE', 'metH', 'lysC', 'thrA', 'metL', 'asd'}

# Fill in gene counts to output df
for index, row in df.iterrows():
    name = row['gene_list']
    gen = row['time']
    if name in emp_set:
        df_output.loc[df_output['Generation']==gen,'EMP']+=1
    if name in ed_set:
        df_output.loc[df_output['Generation']==gen,'ED']+=1
    if name in oxppp_set:
        df_output.loc[df_output['Generation']==gen, 'OxPPP']+=1
    if name in noxppp_set:
        df_output.loc[df_output['Generation']==gen, 'NonOxPPP']+=1
    if name in lowgly_set:
        df_output.loc[df_output['Generation']==gen, 'Lower Glycolysis']+=1
    if name in glctrans_set:
        df_output.loc[df_output['Generation']==gen, 'Glucose Transport']+=1
    if name in lgnopyfk_set:
        df_output.loc[df_output['Generation']==gen, 'Lower Glycolysis no pykF']+=1
    if name in met_set:
        df_output.loc[df_output['Generation']==gen, 'Met Biosynth']+=1
    if name in hser_set:
        df_output.loc[df_output['Generation']==gen, 'Homoserine Biosynth']+=1
    if name in met_hser_set:
        df_output.loc[df_output['Generation']==gen, 'Met + Homoser Biosynth']+=1

df_output.to_csv("output_ara-2bMOTs.csv",index=False)