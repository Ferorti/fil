#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from Bio import AlignIO
import pandas as pd
import numpy as np

# alineamiento de entrada de las secuencias utilizadas para caluclar los scores,
# de aqui se sacan las posiciones que tienen gaps
entrada = "/home/fernando/git/flavivirus-disorder/2019/flavivirus/no_X/hit08/poliproteina_alineado.fasta"
tipo = "fasta"
ali = AlignIO.read(entrada, tipo)

# se lee el archivo de scores, para cada secuencia una linea, con scores separados por comas.
with open("/home/fernando/restaurado/alineamiento/nuevos_scores") as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content]

#convierte a una array.
series = []
for i in range(len(ali)):
    series.append(np.array(content[i].split(",")))

final = []
final2 = []

ID = 0
for virus in ali:
    print virus.id
    secuencia = virus.seq
    print secuencia

    scores_seq = np.array([])
    ls = []
    pos_score = 0
    print len(secuencia)
    for residuo in secuencia:
        #print residuo

        if residuo != "-":
            #ls.append(series[ID][pos_score)

            scores_seq = np.append(scores_seq, series[ID][pos_score])
            #print ID, pos_score,  series[ID][pos_score]
            pos_score = pos_score + 1
        else:
            scores_seq = np.append(scores_seq, -1)
            ls.append(-1)
    ID = ID + 1
    final.append(scores_seq)
#    final2.append(ls)


df = pd.DataFrame(list(map(np.ravel, final)))

df.to_csv("/home/fernando/restaurado/alineamiento/scores_df2", sep='\t', index=False)


