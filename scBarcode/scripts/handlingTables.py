import pandas as pd
import os 
import numpy as np
import argparse
#add arguments parser
parser = argparse.ArgumentParser(description='Handling demultiplexed BAM files')
parser.add_argument('-i', '--input', help='Inpurt you demutplexed Best Score', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)   
args = parser.parse_args()



CellIDs_original = pd.read_csv(args.input , sep="\t" )
ine_dict = {"KOLF2C1day25": "CTL08A", "S20201":"CTL04E", "MIFF1day25":"CTL02A"}
CellID = CellIDs_original.query("Consensus not in ['doublet' , 'LowQuality']")
CellID["Consensus"] = CellID["Consensus"].replace(ine_dict)

for genotype in CellID.Consensus.unique():
    #write index as txt file 
    index = CellID.query("Consensus == @genotype").barcode.tolist()  # Get the list of barcodes for the current genotype
    with open(f"{args.output}/barcodes/{genotype}_barcodes.txt", "w") as f:
        f.write("\n".join(index))  # Write the barcodes to a text file, one per line    