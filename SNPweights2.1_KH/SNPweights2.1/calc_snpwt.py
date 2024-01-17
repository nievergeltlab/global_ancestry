#!/usr/bin/python
from optparse import OptionParser
import numpy
import math
def calculate_snpwt():
	# parse command
    parser = OptionParser()
    parser.add_option("-p", "--par", dest="input_par", metavar="FILE", help="input parameter file path")
    (options, args) = parser.parse_args()
    input_par = options.input_par

    # input parameter
    par_file = open(input_par, "r")
    input = par_file.readlines()
    for line in input:
        if line.split(":")[0]=="geno":
            geno_file = open(line.split()[1],"r")
        if line.split(":")[0]=="snp":
            snp_file = open(line.split()[1],"r")
        if line.split(":")[0]=="ind":
            ind_file = open(line.split()[1],"r")
        if line.split(":")[0]=="evec":
            evec_file = open(line.split()[1],"r")
        if line.split(":")[0]=="eval":
            eval_file = open(line.split()[1],"r")
        if line.split(":")[0]=="snpwtoutput":
            output_file = open(line.split()[1],"w")
        if line.split(":")[0]=="log":
            log_file = open(line.split()[1],"r")

    geno = geno_file.readlines()
    snp  = snp_file.readlines()
    ind  = ind_file.readlines()
    evec = evec_file.readlines() 
    eval = eval_file.readlines()   
    log  = log_file.readlines()

    # input parameters
    eval_list = []
    for i in range(len(eval)):
        eval_list.append(float(eval[i]))# EVAL
    nsnp = float(len(snp)) # nSNP
    nsample = int(len(geno[0])-1) # N_Sample
    iteration = 200 # N.Iter

    # create poplist and npc
    poplist_ind = []
    for row in range(len(ind)):
        poplist_ind.append(ind[row].rstrip().split()[2])
    popset = set(poplist_ind)
    poplist = list(popset)
    npc = len(poplist)-1
    
    # calcualte shrinkage
    gamma1 = nsnp / float(nsample)
    l_pred_a = [0.0]*nsample
    sum_eval = sum(eval_list)
        
    for i in range(iteration):
        number = 0
        for j in range(nsample):
            l_obs = eval_list[j] / sum_eval * (nsnp + sum(l_pred_a)) 
            if (math.pow((1.0+l_obs-gamma1),2)-4.0*l_obs) < 0.0:
                number = j
                break
        if number == 0:
            number = 1
        for k in range(number):
            l_obs = eval_list[k] / sum_eval * (nsnp + sum(l_pred_a)) 
            l_pred = ((1+l_obs - gamma1)+ math.sqrt(math.pow((1.0+l_obs-gamma1),2)-4.0*l_obs))/2.0
            l_pred_a[k] = l_pred
    
    for i in range(npc):            
        l_obs = eval_list[i] / sum_eval * (nsnp + sum(l_pred_a)) 
        l_pred = ((1.0+l_obs-gamma1)+ math.sqrt(math.pow((1.0+l_obs-gamma1),2)-4.0*l_obs))/2.0
        a = math.sqrt((1.0-(gamma1/((l_pred-1.0)**2)))/(1.0+(gamma1/(l_pred-1.0))))
        factor = 1.0/(math.sqrt(l_obs)/math.sqrt(l_pred*(a**2)-(a**2)+1.0))
        print >>output_file, factor,
    print >>output_file, ""

    # calcualte linear transformation parameters
    for i in poplist:
        print >>output_file, i,
    print >>output_file, ""

    popcount = []
    for i in poplist:
        popcount.append(poplist_ind.count(i))
        print >>output_file, poplist_ind.count(i),
    print >>output_file, ""

    popdic = {}
    for i in range(len(poplist)):
        popdic[poplist[i]] = i

    evec_popsum = [[0.0]*npc]*len(poplist)
    for row in range(1,len(evec)):
        evec_tmp = evec[row].rstrip().split()[1:(npc+1)]
        evec_tmp = [float(x) for x in evec_tmp]
        pop_tmp = poplist_ind[(row-1)]
        evec_popsum[popdic[pop_tmp]] = [a + b for (a,b) in zip(evec_tmp,evec_popsum[popdic[pop_tmp]])]
    evec_popavg = []
    for i in poplist:
        evec_popavg.append([a / b for (a,b) in zip(evec_popsum[popdic[i]],[popcount[popdic[i]]]*npc)]+[1.0])

    for i in range(len(evec_popavg)):
        print >>output_file, poplist[i],
        for j in evec_popavg[i][0:-1]:
            print >>output_file, j,
    print >>output_file, ""

    # solve for coeffecients converting PCs to % ancestry 
    # assuming each population group in the reference samples corresponses to 100% ancestry for that particular population
    poppar = []  
    for i in range(len(poplist)):
        b_list = [0.0]*len(poplist)  
        b_list[i] = 1.0
        b = numpy.array(b_list)
        a = numpy.array(evec_popavg)
        x = numpy.linalg.solve(a,b)
        poppar.append(x)

    # solve for coeffecients specific to the reference samples of YRI, CEU, ASI (CHB+CHD) (from HapMap 3), NAT (Reich et al. 2012 Nature)
    # poppar = []  
    # big_b_list = [[1.0,0.0,0.0,0.007206322],[0.0,1.0,0.0,0.060966667],[0.0,0.0,1.0,0.0],[0.0,0.0,0,0.931827]]
    # for i in range(len(poplist)):
    #     b_list = big_b_list[i]
    #     b = numpy.array(b_list)
    #     a = numpy.array(evec_popavg)
    #     x = numpy.linalg.solve(a,b)
    #     poppar.append(x)

    for i in poppar:
        for j in i:
            print >>output_file, j,
    print >>output_file, ""

    # calculate avrage allele count for each SNP
    avg_allele_count = []
    for row in range(len(geno)):
        genotype = geno[row].split()[0]
        snp_id = snp[row].split()[0]
        n = 0
        k = 0
        for col in range(len(genotype)):
            if int(genotype[col]) == 9: continue    
            else:
                n += 1
                k = k + int(genotype[col])
        avg_allele_count.append(float(k)/float(n))
        
    # extract tr(XTX)/(n-1)
    for line in range(len(log)):
        if log[line][0:6] == "trace:":
            trace = float(log[line].split()[1])
    
    # calculate SNP weights U = XVSi
    for row in range(len(geno)): 
        xmean = float(avg_allele_count[row])
        if 2.0 > xmean > 0.0:
            snp_info = snp[row].split()
            print >>output_file, snp_info[0], snp_info[4], snp_info[5], ("%.5f" % (xmean/2.0)),
            genotype = geno[row].rstrip()
            for pc in range(npc):
                S2 = float(evec[0].split()[(pc+1)])
                XV = 0.0
                for col in range(len(genotype)):
                    geno_tmp = float(genotype[col])
                    if geno_tmp == 9:
                        geno_tmp = 0.0
                    XV = XV + float(evec[col+1].split()[pc+1])*((geno_tmp-xmean)/(math.sqrt((xmean/2.0)*(1.0-(xmean/2.0)))))
                W = XV/(S2*trace)
                if pc == (npc-1):
                    print >>output_file, ("%.4e" % W)
                else:
                    print >>output_file, ("%.4e" % W),
        else: continue
                
    geno_file.close()
    snp_file.close()
    ind_file.close()
    evec_file.close()
    eval_file.close()
    output_file.close()
    par_file.close()
       
if __name__ == "__main__":
    calculate_snpwt()