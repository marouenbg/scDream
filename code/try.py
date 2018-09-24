#Script for the single cell dream challenge
import os#change working dir
import pandas as pd

#change working directory
path='/home/marouen/dreamChallenge'
os.chdir(path)

#load data
sc=pd.read_csv('data/singleCellData/dge_normalized.txt',delimiter='\t')
ref=pd.read_csv('data/refDB/bdtnp.txt',delimiter='\t')