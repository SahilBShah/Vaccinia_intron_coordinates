#Package imports
import AGEpy as age
import argparse
import copy
from gtfparse import read_gtf
import numpy as np
import pandas as pd

def main():

	#Command line argument for inputted gtf file
	parser = argparse.ArgumentParser()
	gtf_file = parser.add_argument('gtf_file', type=str, help='Input the gtf file containing genetic coordinates from RNA-seq data. Make sure the gtf extension is included.')
	args = parser.parse_args()

	#Read in gtf file and format necessary columns
	genes_df = read_gtf(args.gtf_file)
	genes_df_frames = genes_df[genes_df["feature"] != "CDS"]
	genes_df_frames["frame"] = "."
	genes_df_cds = genes_df[genes_df["feature"] == "CDS"]
	genes_df = pd.concat([genes_df_frames, genes_df_cds])
	genes_df = genes_df.sort_values(by="seqname")
	genes_df["score"] = "."

	#Filter for only exons
	genes_df_cp = copy.deepcopy(genes_df)
	genes_df_cp = genes_df_cp[genes_df_cp["feature"] == "exon"]

	#Drop unnecessary columns
	genes_df_cp = genes_df_cp.drop(columns=["seqname", "source", "score", "frame", "gene_name", "transcript_id", "tss_id", "p_id"])

	#Dropping duplicate rows
	genes_df_cp = genes_df_cp.drop_duplicates()

	#Obtain all unique genes
	unique_genes = genes_df_cp["gene_id"].unique().tolist()
	unique_genes.sort()

	#Begin intron analysis
	#Identify introns for every gene
	#This may take a few minutes
	df_lst = []
	gene_starts = genes_df[genes_df["feature"] == "start_codon"].reset_index().drop(columns=["index"])
	gene_stops = genes_df[genes_df["feature"] == "stop_codon"].reset_index().drop(columns=["index"])

	for gene in unique_genes:
	    print(gene)
	    one_gene_df = genes_df_cp[genes_df_cp["gene_id"] == gene].sort_values(by="start")
	    
	    if len(one_gene_df) > 1:
	        adj_df = one_gene_df.sort_values(by=["gene_id", "start"]).drop(index=one_gene_df.index[0], axis=0, inplace=False)

	        new_starts = adj_df["start"].to_list()
	        if not gene_stops[gene_stops["gene_id"] == gene].empty and gene_stops[gene_stops["gene_id"] == gene]["start"].iloc[0] > max(one_gene_df["end"]):
	            new_starts.append(gene_stops[gene_stops["gene_id"] == gene]["start"].iloc[0])
	        else:
	            one_gene_df.drop(one_gene_df.tail(1).index,inplace=True)

	        one_gene_df["Next_exon_start"] = new_starts
	        #Remove exons that have the same start and next exon start
	        lowest_idx = one_gene_df.index[one_gene_df["start"] == one_gene_df.iloc[-1]["Next_exon_start"]].tolist()
	        if len(lowest_idx) > 0:
	            for i in lowest_idx:
	                if one_gene_df.at[i, "start"] == one_gene_df.at[i, "Next_exon_start"]:
	                    one_gene_df = one_gene_df.drop([i])
	        
	        if len(one_gene_df) > 0:
	            #If the exon starting positions are the same for the last exons
	            last_next_exon = one_gene_df.iloc[-1]["Next_exon_start"]
	            last_exon_starts_lst = one_gene_df.index[one_gene_df["start"] == one_gene_df.iloc[-1]["start"]].tolist()
	            if last_exon_starts_lst[-1] == one_gene_df.last_valid_index():
	                last_exon_starts_lst.pop()
	            if len(last_exon_starts_lst) > 0:
	                for idx in last_exon_starts_lst:
	                    one_gene_df.at[idx, "Next_exon_start"] = last_next_exon

	            #one_gene_df["Intron_start"] = one_gene_df["Next_exon_start"] - abs(one_gene_df["end"] - one_gene_df["Next_exon_start"]) + 1 
	            one_gene_df["Intron_start"] = one_gene_df["end"] + 1
	            one_gene_df["Intron_stop"] = one_gene_df["Next_exon_start"] - 1
	            one_gene_df = one_gene_df.reset_index().drop(columns=["index"])

	            if not gene_starts[gene_starts["gene_id"] == gene]["start"].empty:
	                if min(one_gene_df["start"]) != gene_starts[gene_starts["gene_id"] == gene]["start"].iloc[0] and gene_starts[gene_starts["gene_id"] == gene]["start"].iloc[0] < min(one_gene_df["start"]):
	                    initial_intron_df = pd.DataFrame({"feature": ["intron"], "start": [np.nan], "end": [np.nan],
	                                                      "strand": [one_gene_df["strand"].iloc[0]],
	                                                      "gene_id": [gene], "Next_exon_start": [min(one_gene_df["start"])], 
	                                                       "Intron_start": [gene_starts[gene_starts["gene_id"] == gene]["end"].iloc[0] + 1],
	                                                      "Intron_stop": [min(one_gene_df["start"]) - 1]})
	                    one_gene_df.append(initial_intron_df)

	            #Deals with alternatively spliced exons
	            grouped_df = one_gene_df.groupby("start").mean().reset_index()
	            same_starts = (one_gene_df["start"] == one_gene_df["Next_exon_start"])
	            same_starts_df = one_gene_df.loc[same_starts, ["start", "Next_exon_start"]]["start"]
	            if len(same_starts_df) > 0:
	                for current_start in list(set(same_starts_df.tolist())):
	                    for i in one_gene_df.index[one_gene_df["start"] == current_start].tolist():
	                        if not grouped_df.iloc[grouped_df.index[grouped_df['start'] == current_start].tolist()[0]].equals(grouped_df.iloc[-1]):
	                            start_replace = grouped_df.iloc[grouped_df.index[grouped_df['start'] == current_start].tolist()[0] + 1]["start"]
	                            one_gene_df.at[i, "Next_exon_start"] = start_replace
	                            one_gene_df.at[i, "Intron_stop"] = one_gene_df.iloc[i]["Next_exon_start"] - 1

	            df_lst.append(one_gene_df)
	    else:
	        if not gene_starts[gene_starts["gene_id"] == gene].empty:
	            if one_gene_df["start"].iloc[0] != gene_starts[gene_starts["gene_id"] == gene]["start"].iloc[0] and gene_starts[gene_starts["gene_id"] == gene]["end"].iloc[0] < one_gene_df["start"].iloc[0]:
	                initial_intron_df = pd.DataFrame({"feature": ["intron"], "start": [np.nan], "end": [np.nan],
	                                              "strand": [one_gene_df["strand"].iloc[0]],
	                                              "gene_id": [gene], "Next_exon_start": [min(one_gene_df["start"])], 
	                                               "Intron_start": [gene_starts[gene_starts["gene_id"] == gene]["end"].iloc[0] + 1],
	                                              "Intron_stop": [min(one_gene_df["start"]) - 1]})
	                df_lst.append(initial_intron_df)
	        if not gene_stops[gene_stops["gene_id"] == gene].empty:
	            if one_gene_df["end"].iloc[0] != gene_stops[gene_stops["gene_id"] == gene]["end"].iloc[0] and gene_stops[gene_stops["gene_id"] == gene]["start"].iloc[0] > one_gene_df["end"].iloc[0]:
	                initial_intron_df = pd.DataFrame({"feature": ["intron"], "start": [np.nan], "end": [np.nan],
	                                                  "strand": [one_gene_df["strand"].iloc[0]],
	                                                  "gene_id": [gene], "Next_exon_start": [gene_stops[gene_stops["gene_id"] == gene]["start"]].iloc[0], 
	                                                   "Intron_start": [max(one_gene_df["start"]) + 1],
	                                                  "Intron_stop": [gene_stops[gene_stops["gene_id"] == gene]["start"].iloc[0] - 1]})
	                df_lst.append(initial_intron_df)
	        
	    
	#Uncomment if introns are needed when end positions are the same######################################################
	#    #Obtains duplicated end values
	#     same_ends_df = one_gene_df.loc[one_gene_df[one_gene_df["gene_id"] == gene].duplicated(subset=["end"], keep=False)]
	#     for unique_end in same_ends_df["end"].unique().tolist():
	#         #Subset df for each same end position
	#         tmp2 = one_gene_df[one_gene_df["end"] == unique_end]
	#         #ID first instance of preceding value
	#         preceding_row = one_gene_df.iloc[[tmp2.index[0] - 1]]
	#         if 0 not in tmp2.index.tolist():
	#             #Drop first because we already have it
	#             tmp2.drop(tmp2.index[0], inplace=True)
	#             #Add preceding row to df
	#             tmp2.append(preceding_row)
	#             #Sort values by start so preceding row is on top
	#             tmp2 = tmp2.sort_values(by="start")
	#             #Get new intron start and end positions
	#             intron_end_lst = [start-1 for start in tmp2["start"].tolist()]
	#             end_length = len(intron_end_lst)
	#             intron_start_lst = [preceding_row["Intron_start"].iloc[0]] * end_length
	#             add_df = pd.DataFrame({"feature": [np.nan]*end_length, "start": [np.nan]*end_length, "end": [np.nan]*end_length,
	#                                                       "strand": [np.nan]*end_length,
	#                                                       "gene_id": [gene]*end_length, "Next_exon_start": [np.nan]*end_length, 
	#                                                        "Intron_start": intron_start_lst,
	#                                                       "Intron_stop": intron_end_lst})
	#             df_lst.append(add_df)
	######################################################################################################################   
	    
	exons_introns_df = pd.concat(df_lst)
	index_names = (exons_introns_df["Intron_start"] < exons_introns_df["Intron_stop"])
	exons_introns_df = exons_introns_df.loc[index_names]

	#Merge dataframes containing all the necessary information to create a gtf file
	to_merge_df = genes_df[["seqname", "source", "feature", "strand", "frame", "gene_id", "gene_name"]]
	to_merge_df = to_merge_df[to_merge_df["feature"] == "exon"].drop(columns="feature")
	to_merge_df = to_merge_df.drop_duplicates()

	#Create dataframe only containing introns and format as a gtf file
	introns_df = pd.DataFrame({"start": exons_introns_df["Intron_start"].tolist(), "end": exons_introns_df["Intron_stop"].tolist(), 
	                          "gene_id": exons_introns_df["gene_id"].tolist()})
	introns_df["feature"] = "intron"
	introns_df["score"] = "."
	introns_df["transcript_id"] = np.nan
	introns_df["tss_id"] = np.nan
	introns_df["p_id"] = np.nan
	introns_df = introns_df.merge(to_merge_df, on="gene_id", how="left")
	introns_df = introns_df[["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id",
	                       "gene_name", "transcript_id", "tss_id", "p_id"]]
	introns_df["frame"] = "."
	introns_df.replace(np.nan, "", inplace=True)

	#Create the final dataframe
	final_df = pd.concat([genes_df, introns_df])

	#Output final dataframe as a gtf file
	age.writeGTF(final_df,"./genes_w_introns.gtf")

if __name__ == '__main__':
    main()








