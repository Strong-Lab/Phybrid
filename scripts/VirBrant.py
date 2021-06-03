#!/usr/bin/env python
import argparse
import os
import sys
import glob
import pandas as pd
import numpy as np
import subprocess
import joblib
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Get Kmers and Format into Dictionary')
parser.add_argument('-i','--input', type=str, help='Fasta File')
parser.add_argument('-o','--output', type=str, help="Output_Directory", default="VirBrant_Output")
parser.add_argument('-m','--model', type=str, help="Model files")
args = parser.parse_args()


command_dir = "mkdir -p " + args.output
os.system(command_dir)

command = "awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' " + args.input + " > " + args.output + "/Sizes.txt"
os.system(command)

def filter_fasta(fasta_file, list_of_headers, outputs):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open(outputs, "w") as out_file:
        for fasta in fasta_sequences:
            if fasta.description in list_of_headers:
                SeqIO.write(fasta, out_file, "fasta")
            else:
                continue
                
def filter_fasta_2(fasta_file, list_of_headers, outputs):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open(outputs, "w") as out_file:
        for fasta in fasta_sequences:
            if fasta.id in list_of_headers:
                SeqIO.write(fasta, out_file, "fasta")
            else:
                continue


### Fasta Split
## Split Fasta
def fasta_split(input_fasta, size_file):

    def dictionary_to_list_header(dictionary):
        my_list = []
        for item in dictionary.keys():
            my_list.append(item.split(">")[1])
        return(my_list)


    indecies = []
    lengths = []
    with open(size_file, "r") as my_file:
        for index, line in enumerate(my_file.readlines()):
            line = line.strip()
            if index%2==0:
                indecies.append(line)
            else:
                lengths.append(int(line))
    my_file.close()

    dictionary = dict(zip(indecies, lengths))

    dictionary_1K = {k: v for k, v in dictionary.items() if v >= 750 and v < 2000}
    dictionary_3K = {k: v for k, v in dictionary.items() if v >= 2000 and v < 4000}
    dictionary_5K = {k: v for k, v in dictionary.items() if v >= 4000 and v < 7500}
    dictionary_10K = {k: v for k, v in dictionary.items() if v >= 8000}

    ## 1K
    if len(dictionary_1K)>0:
        outp = dictionary_to_list_header(dictionary_1K)
        output_loc = args.output+"/temp_1k.fasta"
        filter_fasta(input_fasta, outp, output_loc)
    
    if len(dictionary_3K)>0:
        outp = dictionary_to_list_header(dictionary_3K)
        output_loc = args.output+"/temp_3k.fasta"
        filter_fasta(input_fasta, outp, output_loc)
    
    ## 5K
    if len(dictionary_5K)>0:
        outp = dictionary_to_list_header(dictionary_5K)
        output_loc = args.output+"/temp_5k.fasta"
        filter_fasta(input_fasta, outp, output_loc)
    
    ## 10K
    if len(dictionary_10K)>0:
        outp = dictionary_to_list_header(dictionary_10K)
        output_loc = args.output+"/temp_10k.fasta"
        filter_fasta(input_fasta, outp, output_loc)
    

sizes = args.output+"/Sizes.txt"
fasta_split(args.input, sizes)

x = pd.read_csv("Full_Features.txt", header=None)
columns_to_keep = list(x[0])

##### Kmerizer
def kmerizer(fastaFile, k, output_size):
    kmerCmd = 'scripts/kmer-counter-master/kmer-counter --fasta --k=%d %s' % (k, fastaFile)

    try:
        output = subprocess.check_output(kmerCmd, shell=True)
        output = output.decode()
        result = {}
        for line in output.splitlines():
            header, counts = line.strip().split('\t')
            header = header[1:]
            kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in counts.split(' ')])
            result[header] = kmers

        diction = pd.DataFrame.from_dict(result, orient="index")
        diction.replace(np.nan, 0, inplace=True)
        mytrans = str.maketrans('ATCG', 'TAGC')
        ## Merge Compliments
        diction.columns = np.sort([diction.columns, [x.translate(mytrans) for x in diction.columns]], axis=0)[0, :]
        diction = diction.groupby(level=0, axis=1).sum()
        mytrans = str.maketrans('ATCG', 'TAGC')
        diction.columns = np.sort([diction.columns, [x.translate(mytrans)[::-1] for x in diction.columns]], axis=0)[0, :]
        diction = diction.groupby(level=0, axis=1).sum()
        diction.columns = np.sort([diction.columns, [x[::-1] for x in diction.columns]], axis=0)[0, :]
        diction = diction.groupby(level=0, axis=1).sum()
        diction = diction.reindex(columns_to_keep, axis=1)
        
        
        
        ## Normalize by Length
        diction = diction.div(diction.sum(axis=1), axis=0)
        diction = diction.round(7)
    
        name = output_size + ".csv"
        diction.to_csv(name)
    
    except subprocess.CalledProcessError as error:
        sys.stderr.write("%s\n" % (str(error)))


#### Machine Learning Predictions
def machine_learnizer(model, input_file, input_fasta, proteins=0):
    df = pd.read_csv(input_file)
    df['Unnamed: 32774'] = 0
    ## Move ID to Index
    try:
        df = df.set_index('Unnamed: 0')
    except:
        df = df.set_index("IDs")

    ## Make sure kmer columns match models
    z = pd.read_csv("scripts/Models/Headers.txt", header=None)
    q = z.transpose()
    cols = set(q[0])
    colums_to_add = set(cols) - set(df.columns)
    if len(colums_to_add) > 0:
        for item in colums_to_add:
            df[item] = 0.0
            
    
    ## Load model and Predict on dataset
    c = joblib.load(model)
    names = c.get_booster().feature_names
    q = c.predict(df[names])
    
    output_df = pd.DataFrame()
    output_df['Contigs'] = df.index
    output_df["Predicted_Values"] = list(q)
    output_df = output_df[output_df["Predicted_Values"]=="Phage"]
    phage_output = list(output_df['Contigs'])
    
    
    output_file = input_file.split(".")[0]
    output_file = output_file + "_viral.fasta"
    if proteins==0:
        filter_fasta(input_fasta, phage_output, output_file)
    else:
        filter_fasta_2(input_fasta, phage_output, output_file)


def Operon_Distribution(list_of_strains):
    """
    Input: Pandas Strand Column
    Output: Number of Operons, mean size of operons, stdev of operons, median of operons
    Usage:
    
    """
    gene_strand_list = []
    strands = list(list_of_strains)
    counts = 1
    if len(strands)==1:
        return(len(strands), len(strands), len(strands), len(strands))
    for i in range(0, len(strands)-1):
        if strands[i]==strands[i+1]:
            counts+=1
        
        ## Deal with End Case
            if (i+2)==len(strands):
                gene_strand_list.append(counts)
        else:
            gene_strand_list.append(counts)
            counts = 1
    return(len(gene_strand_list), max(gene_strand_list), np.mean(gene_strand_list), np.median(gene_strand_list))


def prodigal_parser(fasta_size, protein_size):
    f = open(fasta_size, "r")
    lengths = f.read().split(">")[1:]
    f.close()
    df = pd.DataFrame(lengths)
    df['ID'] = [x.split("\n")[0] for x in df[0]]
    df['Size'] = [x.split("\n")[1] for x in df[0]]
    df_lengths = df[["ID", "Size"]]
    
    p = open(protein_size, "r")
    lengths_p = p.read().split(">")[1:]
    p.close()
    df_p = pd.DataFrame(lengths_p)
    df_p['ID'] = [x.split("\n")[0] for x in df_p[0]]
    df_p['Size'] = [x.split("\n")[1] for x in df_p[0]]
    df = df_p[["ID", "Size"]]
    
    Gene_Count = len(df)
    ##########################################################################
    ## Hard coded will need to change | Underscores specific to current output
    df['Genome'] = ["_".join(x.split('#')[0].rstrip().split("_")[0:-1]) for x in df["ID"]]
    df['Strand'] = [pd.to_numeric(x.split('#')[3].rstrip().lstrip()) for x in df["ID"]]
    df['GC'] = [pd.to_numeric(x.split('#')[4].split(';')[5].split("=")[1]) for x in df["ID"]]
    df['Start'] = [pd.to_numeric(x.split('#')[1].rstrip().lstrip()) for x in df["ID"]]
    df['Stop'] = [pd.to_numeric(x.split('#')[2].rstrip().lstrip()) for x in df["ID"]]
    df.drop("ID", axis=1, inplace=True)
    Groups_df = [v for k, v in df.groupby('Genome', as_index=False)]
    ###########################################################################
    ## This chunk creates a dictionary containing genome strand absolute values
    genome_strands = {}
    for file in range(len(Groups_df)):
        df2 = pd.DataFrame(Groups_df[file])
        Gene_Count = len(df2)

        ## Get GC Percentage
        GC = df2["GC"].median()
    
        AA_Length = df2["Size"].median()
    
        ## Get Overlap Percentage
        Overlap = 0
        if len(df2)==1:
            Percent_Overlap=0
        else:
            for i in range(len(df2)-1):
                if df2.iloc[i]['Strand'] != df2.iloc[i+1]['Strand']:
                    next
                else:
                    if df2.iloc[i]['Stop'] > df2.iloc[i+1]['Start']:
                        Overlap += 1
            Percent_Overlap = Overlap/len(df2)
    
    
        ids = list(df2["Genome"])[0]
        length = df_lengths[df_lengths['ID'].str.match(ids)]

        if len(length)==0:
        #print("Empty: ")
            continue
        lengths = int(length.iloc[0][1])
        Gene_Density = float(Gene_Count)/(lengths/1000)
        num_operons, max_operon, mean_operon_length, median_operon_length = Operon_Distribution(df2["Strand"])
    
        genome_strands[ids]=[Gene_Density, median_operon_length, Percent_Overlap, AA_Length]
    
    df_features = pd.DataFrame(genome_strands)
    df_features = df_features.transpose()
    df_features.columns = ["Gene_Density", "Median Operon Length", "Percent_Overlap", "AA_Length"]
    Groups_df=1
    return(df_features)
        
def proteinizer(input_file, output_file, fasta_lengths, kmer_csv):
    command = "prodigal -i " + input_file + " -a " + output_file + ".faa -p meta"
    os.system(command)
    proteins = output_file + ".faa"
    command = "awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' " + proteins+ " > " + output_file + ".txt"
    os.system(command)
    protein_lengths = output_file + ".txt"
    q = prodigal_parser(fasta_lengths, protein_lengths)
    kmer_df = pd.read_csv(kmer_csv)
    kmer_df = kmer_df.reset_index()
    kmer_df['IDs'] = [x.split(" ")[0] for x in kmer_df["Unnamed: 0"]]
    kmer_df = kmer_df.set_index("IDs")
    kmer_df = kmer_df.drop(["Unnamed: 0","index"], axis=1)
    kmer_df = pd.merge(kmer_df, q, left_index=True, right_index=True)
    kmer_df.to_csv(kmer_csv)


directory_path = os.getcwd()


outf = args.output+"/temp_1k.fasta"
outs = args.output+"/temp_1k"
outp = args.output+"/temp_1k.csv"
outt = args.output+"/temp_1k.txt"
kmerizer(outf, 8, outs)
proteinizer(outf, outs, sizes, outp)
model_path_1k = directory+"/Models/Features_1K.cls"
machine_learnizer(model_path_1k,outp, outf, proteins=1)
rm_command = "rm " + outf + " " + outp + " "+ outs + ".faa " + outt
os.system(rm_command)


outf = args.output+"/temp_3k.fasta"
outs = args.output+"/temp_3k"
outp = args.output+"/temp_3k.csv"
outt = args.output+"/temp_3k.txt"
kmerizer(outf, 8, outs)
proteinizer(outf, outs, sizes, outp)
model_path_3k = directory+"/Models/Features_3K.cls"
machine_learnizer(model_path_3k,outp, outf, proteins=1)
rm_command = "rm " + outf + " " + outp + " "+ outs + ".faa " + outt
os.system(rm_command)

outf = args.output+"/temp_5k.fasta"
outs = args.output+"/temp_5k"
outp = args.output+"/temp_5k.csv"
outt = args.output+"/temp_5k.txt"
kmerizer(outf, 8, outs)
proteinizer(outf, outs, sizes, outp)
model_path_5k = directory+"/Models/Features_5K.cls"
machine_learnizer(model_path_5k,outp, outf, proteins=1)
rm_command = "rm " + outf + " " + outp + " " + outs + ".faa " + outt
os.system(rm_command)

outf = args.output+"/temp_10k.fasta"
outs = args.output+"/temp_10k"
outp = args.output+"/temp_10k.csv"
outt = args.output+"/temp_10k.txt"
kmerizer(outf, 8, outs)
proteinizer(outf, outs, sizes, outp)
model_path_10k = directory+"/Models/Features_10K.cls"
machine_learnizer(model_path_10k,outp, outf, proteins=1)

rm_command = "rm " + outf + " " + outp + " " + outs + ".faa " + outt
os.system(rm_command)


command = "cat " + args.output+"/*_viral.fasta > " + args.output + "/VirBrant_Viral_Predictions.fasta"
os.system(command)
os.system("rm " + args.output+ "/*_viral.fasta")
os.system("grep '>' "+args.output+"/VirBrant_Viral_Predictions.fasta > "+args.output+"/VirBrant_Headers.txt")
