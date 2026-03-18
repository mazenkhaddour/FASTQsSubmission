import pandas as pd
import os 
import numpy as np
import argparse
#add arguments parser
parser = argparse.ArgumentParser(description='Handling demultiplexed BAM files')
parser.add_argument('-i', '--input', help='Inpurt you demutplexed Best Score', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)   
args = parser.parse_args()


CellIDs = pd.read_csv(args.input , sep="\t" )
Singlet = CellIDs.query("DropletType == 'Singlet'")  # Filter the CellIDs DataFrame to include only rows where DropletType is 'Singlet'
Singlet_barcodes = Singlet.loc[:,["barcode","FirstID"]]
#creat folder to store 
if os.path.exists(args.output+"/barcodes") == False:
    os.makedirs(args.output+"/barcodes")
for genotype in Singlet_barcodes.FirstID.unique():
    #write index as txt file 
    index = Singlet_barcodes.query("FirstID == @genotype").barcode.tolist()  # Get the list of barcodes for the current genotype
    with open(f"{args.output}/barcodes/{genotype}_barcodes.txt", "w") as f:
        f.write("\n".join(index))  # Write the barcodes to a text file, one per line    