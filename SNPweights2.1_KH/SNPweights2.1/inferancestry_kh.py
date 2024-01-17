#!/usr/bin/python
from heapq import merge
from optparse import OptionParser
import numpy as np
import pandas as pd

def infer_anc():
    # parse command
    parser = OptionParser()
    parser.add_option("-p", "--par", dest="par", metavar="FILE", help="input parameter file path")
    (options, args) = parser.parse_args()
    input_par = options.par

    # input paramenter and files
    input_file = open(input_par, "r")
    input = input_file.readlines()
    for line in input:
        if line.split(":")[0]=="geno":
            geno_file = line.split()[1]
        if line.split(":")[0]=="snp":
            snp_file = line.split()[1]
        if line.split(":")[0]=="ind":
            ind_file = line.split()[1]
        if line.split(":")[0]=="snpwt":
            snpwt_file = line.split()[1]
        if line.split(":")[0]=="predpcoutput":
            output_file = line.split()[1]

    # par_file = open("./inferancestry.info","r")
    # wt_file = open("./snpwt."+pop,"r")
    geno = np.genfromtxt(geno_file, delimiter=1)
    geno[np.where(geno==9)]=np.nan

    snp  = np.loadtxt(snp_file, dtype = str)
    ind  = np.loadtxt(ind_file, dtype = str)
    weight = np.loadtxt(snpwt_file, dtype = str, skiprows=5)

    snpwt_file_txt = open(snpwt_file,"r")
    weight_txt = snpwt_file_txt.readlines()

    # assign number of PC to predict
    npc = len(weight_txt[0].rstrip().split())
    npop = npc+1
    shrinkage = [float(x) for x in weight_txt[0].rstrip().split()]
    coef_pop = [float(x) for x in weight_txt[4].rstrip().split()]
    wt_number = float(len(weight_txt)-5)

    # create a list of SNP mismatches
    snp_row_index = pd.DataFrame(snp[:,[0]]).reset_index()
    snp_row_index = snp_row_index.set_axis(['Index', 'SNP'], axis = 1)

    weight_alleles = pd.DataFrame(weight[:,0:3],columns=['SNP','Ref_A1','Ref_A2'])
    allofus_alleles = pd.DataFrame(snp[:,[0,4,5]],columns=['SNP','AOU_A1','AOU_A2'])

    merged_alleles = allofus_alleles.merge(weight_alleles, on = 'SNP', how = 'inner')
    merged_alleles['Ref_A1_switched'] = np.where(merged_alleles['Ref_A1'] == 'A', 'T',
                                                 np.where(merged_alleles['Ref_A1'] == 'T', 'A',
                                                          np.where(merged_alleles["Ref_A1"] == 'G', 'C',
                                                                   np.where(merged_alleles["Ref_A1"] == 'C', 'G', ''))))
    merged_alleles['Ref_A2_switched'] = np.where(merged_alleles['Ref_A2'] == 'A', 'T',
                                                 np.where(merged_alleles['Ref_A2'] == 'T', 'A',
                                                          np.where(merged_alleles["Ref_A2"] == 'G', 'C',
                                                                   np.where(merged_alleles["Ref_A2"] == 'C', 'G', ''))))
    
    merged_alleles['category'] = np.where((merged_alleles['Ref_A1'] == merged_alleles['AOU_A1']) & (merged_alleles['Ref_A2'] == merged_alleles['AOU_A2']), 'Matched',
                                          np.where((merged_alleles['Ref_A1_switched'] == merged_alleles['AOU_A1']) & (merged_alleles['Ref_A2_switched'] == merged_alleles['AOU_A2']), 'Matched',
                                                   np.where((merged_alleles['Ref_A2'] == merged_alleles['AOU_A1']) & (merged_alleles['Ref_A1'] == merged_alleles['AOU_A2']), 'Swapped',
                                                            np.where((merged_alleles['Ref_A2_switched'] == merged_alleles['AOU_A1']) & (merged_alleles['Ref_A1_switched'] == merged_alleles['AOU_A2']), 'Swapped','Removed'))))

    # make a list of all swaps and removals
    swap_list = merged_alleles[merged_alleles['category'] == 'Swapped']['SNP']
    removal_list = merged_alleles[merged_alleles['category'] == 'Removed']['SNP']

    # handle the swaps first
    swap_rows = snp_row_index[snp_row_index['SNP'].isin(swap_list)]['Index']
    geno[swap_rows,:] = 2 - geno[swap_rows,:]

    # handle the row removals next
    snp = pd.DataFrame(snp)
    snp = snp.set_axis(['SNP', 'CHROM', 'VAL', 'OTHER', 'A1', 'A2'], axis = 1)
    snp = snp[~snp['SNP'].isin(removal_list)]

    weight = pd.DataFrame(weight)
    weight = weight.set_axis(['SNP', 'A1', 'A2', 'P', 'POP1', 'POP2', 'POP3', 'POP4', 'POP5'], axis = 1)
    weight = weight[~weight['SNP'].isin(removal_list)]

    removal_rows = snp_row_index[snp_row_index['SNP'].isin(removal_list)]['Index']
    geno = np.delete(geno, removal_rows, axis = 0)

    snp_row_index = snp_row_index[~snp_row_index['SNP'].isin(removal_list)]

    # generate X, P, and W
    weight_filtered = snp_row_index.merge(weight, on = 'SNP', how = 'left')
    P = np.array(weight_filtered['P'].astype(float))
    W_full = np.array(weight_filtered[['POP1','POP2','POP3','POP4','POP5']].astype(float))

    WX_matrix = (geno.T - P * 2) / np.sqrt(P * (1 - P))

    ind = pd.DataFrame(ind)
    ind = ind.set_axis(['ID', 'Sex', 'Case_Control'], axis = 1)

    pc_colnames = []

    for pc in range(npc):
        W = W_full[:, pc]
        
        WX = WX_matrix * W
        n = np.sum(~np.isnan(WX), axis = 1)

        WX_sum = np.nansum(WX, axis = 1)

        pred_pc_adj = WX_sum * (wt_number / n) / shrinkage[pc]
        ind['PC' + str(pc)] = pred_pc_adj

        pc_colnames.append('PC' + str(pc))

    pop_colnames = []

    for pop in range(npop):
        weights = coef_pop[pop*npop:(pop+1) * npop - 1]
        colname = 'POP' + str(pop)
        
        pop_colnames.append(colname)

        ind[colname] = np.sum(np.array(ind)[:,3:3+npc] * weights, axis = 1) + coef_pop[(pop+1) * npop - 1]

        ind.loc[ind[colname] < 0, colname] = 0
        ind.loc[ind[colname] > 100, colname] = 100
    
    ind['sum_denom'] = ind[pop_colnames].sum(axis = 1)

    for pop in range(npop): 
        ind.loc[ind['sum_denom'] > 1, 'POP' + str(pop)] = ind['POP' + str(pop)] / ind['sum_denom']

    ind['n'] = n

    out_colnames  = ['ID', 'Case_Control', 'n']
    out_colnames.extend(pc_colnames)
    out_colnames.extend(pop_colnames)
    
    ind.head()
    
    out = ind[out_colnames]
    out.to_csv(output_file, sep=' ', header = False, index = False)

    snpwt_file_txt.close()

if __name__ == "__main__":
    infer_anc()
