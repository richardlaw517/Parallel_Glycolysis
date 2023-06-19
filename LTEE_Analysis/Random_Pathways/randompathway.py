import pandas as pd
import csv
import random
import numpy as np

# Import CSV of all genes
with open('allgenes.csv', newline='') as f:
    reader = csv.reader(f)
    genes = list(reader)[1:]
    genes = np.array(genes).flatten() # 1D array of genes

# Create random sets of genes
randomPathways = {}
for x in range (1,1001):
    dictkey = 'R' + str(x)
    randomPathways[dictkey] =  set(np.array(random.choices(genes, k=10)))

#define columns
columns = ["Pathway","Ara+1A","Ara+1B","Ara+2A","Ara+2B","Ara+3A","Ara+3B",
            "Ara+4A","Ara+4B","Ara+5A","Ara+5B","Ara+6A","Ara+6B","Ara-1A",
            "Ara-1B","Ara-2A","Ara-2B","Ara-3A","Ara-3B","Ara-4A","Ara-4B",
            "Ara-5A","Ara-5B","Ara-6A","Ara-6B"]

def counter(gen):
    input = 'all_strain_' + gen + '.csv'

    # Import CSV of genes all mutations from each clone/strain in a single generation
    df = pd.read_csv (input)
    
    # Create blank output df
    df_output = pd.DataFrame(0, index = np.arange(0,1000), columns = columns)
    df_output["Pathway"] = randomPathways.keys()
    df_output = df_output.set_index("Pathway")
    #df_output = pd.read_csv ('gene_list_random.csv')

    for index, row in df.iterrows():
        name, strain = row['gene_list'], row['population']
        for key in randomPathways.keys():
            if name in randomPathways[key]:
                df_output.loc[key][strain] += 1

    output = 'output_randomMOT_' + gen + '.csv'
    df_output.to_csv(output, index=False)
    print('done'+gen)

gens = ['500', '1000', '1500', '2k', '5k', '10k', '15k', '20k', '30k', '40k', '50k']
for gen in gens:
    counter(gen)