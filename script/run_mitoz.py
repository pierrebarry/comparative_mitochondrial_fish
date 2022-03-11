from Bio import SeqIO
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import allel

arg=sys.argv

fasta_sequences = SeqIO.parse(open(arg[2]),'fasta')

pop_sequence = []
seq_sequence = []
name_sequence = []
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	name_sequence += [name]
	seq_sequence += [sequence]
	pop_sequence += [name[5:7]]


gene=pd.DataFrame(
{
	"NAME":name_sequence,
	"POP":pop_sequence,
	"SEQUENCE":seq_sequence
}
)


## Get length
print("Stats")

length=len(gene["SEQUENCE"][0])

biallelic = 0
polyallelic = 0
gaps_all_individuals = 0
gaps = 0

for i in tqdm(range(length)):
	locus = []
	for j in range(gene.shape[0]):
		locus += [gene["SEQUENCE"][j][i]]
	if '-' in locus:
		gaps = gaps + 1
	if len([locus[jj] for jj in range(len(locus)) if locus[jj] == "-"]) == len(locus):
		gaps_all_individuals = gaps_all_individuals + 1
	locus_without_gap=[locus[jj] for jj in range(len(locus)) if locus[jj] != "-"]
	if len(list(np.unique(locus_without_gap))) == 2:
		biallelic = biallelic +1
	elif len(list(np.unique(locus_without_gap))) > 2:
		polyallelic = polyallelic + 1

gc_content_indiv = []
for indiv in tqdm(range(gene.shape[0])):
	tmp = gene["SEQUENCE"][indiv]
	num_gc=len([tmp[jj] for jj in range(len(tmp)) if tmp[jj] =="g" or tmp[jj] == "c"])
	gc_content_indiv += [num_gc/len([tmp[jj] for jj in range(len(tmp)) if tmp[jj] != "-"])]

gc_content = np.mean(gc_content_indiv)	

print("pi-dxy")

num_pi_all = 0
denom_pi_all = 0
num_pi_li = 0
denom_pi_li = 0
num_pi_mu = 0
denom_pi_mu = 0
num_pi_fa = 0
denom_pi_fa = 0
num_pi_ga = 0
denom_pi_ga = 0

num_dxy_li_mu = 0
denom_dxy_li_mu = 0
num_dxy_li_fa = 0
denom_dxy_li_fa = 0
num_dxy_li_ga = 0
denom_dxy_li_ga = 0
num_dxy_mu_fa = 0
denom_dxy_mu_fa = 0
num_dxy_mu_ga = 0
denom_dxy_mu_ga = 0
num_dxy_fa_ga = 0
denom_dxy_fa_ga = 0

for i in tqdm(range(length)):
	locus = []
	for j in range(gene.shape[0]):
		locus += [gene["SEQUENCE"][j][i]]
	
	pop=list(gene["POP"])
	pi_dxy=pd.DataFrame(
	{
		"POP":pop,
		"LOCUS":locus
	}
	)
	pi_dxy = pi_dxy[pi_dxy.LOCUS != "-"]
	if pi_dxy.shape[0] > 1:
		# pi - all
		denom_pi_all = denom_pi_all + (pi_dxy.shape[0]*(pi_dxy.shape[0]-1))/2
		tmp = list(pi_dxy["LOCUS"])
		for first in range(0,len(tmp)-1):
			for second in range(first+1,len(tmp)):
				if tmp[first]!=tmp[second]:
					num_pi_all = num_pi_all + 1
		
		# pi - li
		tmp = list(pi_dxy[pi_dxy.POP == "Li"]["LOCUS"])
		if len(tmp) > 1:
			denom_pi_li = denom_pi_li + (len(tmp)*(len(tmp)-1))/2
			for first in range(0,len(tmp)-1):
				for second in range(first+1,len(tmp)):
					if tmp[first]!=tmp[second]:
						num_pi_li = num_pi_li + 1
		
		if str(arg[1])!="Afall":
			# pi - mu
			tmp = list(pi_dxy[pi_dxy.POP == "Mu"]["LOCUS"])
			if len(tmp) > 1:
				denom_pi_mu = denom_pi_mu + (len(tmp)*(len(tmp)-1))/2
				for first in range(0,len(tmp)-1):
					for second in range(first+1,len(tmp)):
						if tmp[first]!=tmp[second]:
							num_pi_mu = num_pi_mu + 1
		
		# pi - fa
		tmp = list(pi_dxy[pi_dxy.POP == "Fa"]["LOCUS"])
		if len(tmp) > 1:
			denom_pi_fa = denom_pi_fa + (len(tmp)*(len(tmp)-1))/2
			for first in range(0,len(tmp)-1):
				for second in range(first+1,len(tmp)):
					if tmp[first]!=tmp[second]:
						num_pi_fa = num_pi_fa + 1
		
		# pi - ga
		tmp = list(pi_dxy[pi_dxy.POP == "Ga"]["LOCUS"])
		if len(tmp) > 1:
			denom_pi_ga = denom_pi_ga + (len(tmp)*(len(tmp)-1))/2
			for first in range(0,len(tmp)-1):
				for second in range(first+1,len(tmp)):
					if tmp[first]!=tmp[second]:
						num_pi_ga = num_pi_ga + 1
		
		if str(arg[1])!="Afall":
			# dxy - li - mu
			tmp1 = list(pi_dxy[pi_dxy.POP == "Li"]["LOCUS"])
			tmp2 = list(pi_dxy[pi_dxy.POP == "Mu"]["LOCUS"])
			if len(tmp1) >= 1 and len(tmp2) >= 1:
				for first in range(0,len(tmp1)):
					for second in range(0,len(tmp2)):
						if tmp1[first]!=tmp2[second]:
							num_dxy_li_mu = num_dxy_li_mu + 1
							denom_dxy_li_mu = denom_dxy_li_mu + 1
						else:
							denom_dxy_li_mu = denom_dxy_li_mu + 1
		
		# dxy - li - fa 
		tmp1 = list(pi_dxy[pi_dxy.POP == "Li"]["LOCUS"])
		tmp2 = list(pi_dxy[pi_dxy.POP == "Fa"]["LOCUS"])
		if len(tmp1) >= 1 and len(tmp2) >= 1:
			for first in range(0,len(tmp1)):
				for second in range(0,len(tmp2)):
					if tmp1[first]!=tmp2[second]:
						num_dxy_li_fa = num_dxy_li_fa + 1
						denom_dxy_li_fa = denom_dxy_li_fa + 1
					else:
						denom_dxy_li_fa = denom_dxy_li_fa + 1
		
		# dxy - li - ga
		tmp1 = list(pi_dxy[pi_dxy.POP == "Li"]["LOCUS"])
		tmp2 = list(pi_dxy[pi_dxy.POP == "Ga"]["LOCUS"])
		if len(tmp1) >= 1 and len(tmp2) >= 1:
			for first in range(0,len(tmp1)):
				for second in range(0,len(tmp2)):
					if tmp1[first]!=tmp2[second]:
						num_dxy_li_ga = num_dxy_li_ga + 1
						denom_dxy_li_ga = denom_dxy_li_ga + 1
					else:
						denom_dxy_li_ga = denom_dxy_li_ga + 1
		
		if str(arg[1])!="Afall":
			# dxy - mu - fa
			tmp1 = list(pi_dxy[pi_dxy.POP == "Mu"]["LOCUS"])
			tmp2 = list(pi_dxy[pi_dxy.POP == "Fa"]["LOCUS"])
			if len(tmp1) >= 1 and len(tmp2) >= 1:
				for first in range(0,len(tmp1)):
					for second in range(0,len(tmp2)):
						if tmp1[first]!=tmp2[second]:
							num_dxy_mu_fa = num_dxy_mu_fa + 1
							denom_dxy_mu_fa = denom_dxy_mu_fa + 1
						else:
							denom_dxy_mu_fa = denom_dxy_mu_fa + 1
		
		if str(arg[1])!="Afall":
			# dxy - mu - ga
			tmp1 = list(pi_dxy[pi_dxy.POP == "Mu"]["LOCUS"])
			tmp2 = list(pi_dxy[pi_dxy.POP == "Ga"]["LOCUS"])
			if len(tmp1) >= 1 and len(tmp2) >= 1:
				for first in range(0,len(tmp1)):
					for second in range(0,len(tmp2)):
						if tmp1[first]!=tmp2[second]:
							num_dxy_mu_ga = num_dxy_mu_ga + 1
							denom_dxy_mu_ga = denom_dxy_mu_ga + 1
						else:
							denom_dxy_mu_ga = denom_dxy_mu_ga + 1
					
		# dxy - fa - ga
		tmp1 = list(pi_dxy[pi_dxy.POP == "Fa"]["LOCUS"])
		tmp2 = list(pi_dxy[pi_dxy.POP == "Ga"]["LOCUS"])
		if len(tmp1) >= 1 and len(tmp2) >= 1:
			for first in range(0,len(tmp1)):
				for second in range(0,len(tmp2)):
					if tmp1[first]!=tmp2[second]:
						num_dxy_fa_ga = num_dxy_fa_ga + 1
						denom_dxy_fa_ga = denom_dxy_fa_ga + 1
					else:
						denom_dxy_fa_ga = denom_dxy_fa_ga + 1

pi_all = num_pi_all / denom_pi_all
pi_li = num_pi_li / denom_pi_li
if str(arg[1])!="Afall":
	pi_mu = num_pi_mu / denom_pi_mu
pi_fa = num_pi_fa / denom_pi_fa
pi_ga = num_pi_ga / denom_pi_ga

if str(arg[1])!="Afall":
	dxy_li_mu = num_dxy_li_mu / denom_dxy_li_mu
dxy_li_fa = num_dxy_li_fa / denom_dxy_li_fa
dxy_li_ga = num_dxy_li_ga / denom_dxy_li_ga
if str(arg[1])!="Afall":
	dxy_mu_fa = num_dxy_mu_fa / denom_dxy_mu_fa 
	dxy_mu_ga = num_dxy_mu_ga / denom_dxy_mu_ga
dxy_fa_ga = num_dxy_fa_ga / denom_dxy_fa_ga

if str(arg[1])!="Afall":
	da_li_mu = dxy_li_mu - ((pi_li + pi_mu)/2)
da_li_fa = dxy_li_fa - ((pi_li + pi_fa)/2)
da_li_ga = dxy_li_ga - ((pi_li + pi_ga)/2)
if str(arg[1])!="Afall":
	da_mu_fa = dxy_mu_fa - ((pi_mu + pi_fa)/2)
	da_mu_ga = dxy_mu_ga - ((pi_mu + pi_ga)/2)
da_fa_ga = dxy_fa_ga - ((pi_fa + pi_fa)/2)

## Create VCF
print("VCF")
with open("/home/labosea1/MitoZ/ANALYSIS/"+str(arg[1])+"/VCF/vcf_"+str(arg[1])+".vcf", "a") as file_object:
    file_object.write("##fileformat=VCF \n")
    file_object.write("##contig=<ID=MT,length="+str(length)+",species='"+str(arg[1])+">' \n")
    file_object.write("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'> \n")
    file_object.write("#")
    line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
    file_object.write(line)

line_name = ""
for i in range(len(list(gene["NAME"]))):
	line_name+="{}\t".format(gene["NAME"][i])

with open("/home/labosea1/MitoZ/ANALYSIS/"+str(arg[1])+"/VCF/vcf_"+str(arg[1])+".vcf", "a") as file_object:
	file_object.write(line_name)
	file_object.write("\n")

for i in tqdm(range(length)):
	locus = []
	for j in range(gene.shape[0]):
		locus += [gene["SEQUENCE"][j][i]]
	
	locus_without_gap=[locus[jj] for jj in range(len(locus)) if locus[jj] != "-"]
	locus_without_gap=list(np.unique(locus_without_gap))
	if len(locus_without_gap) == 0:
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("MT",str(i),".",".",".",str(10000),"PASS",".","GT")
	elif len(locus_without_gap) == 1:
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("MT",str(i),".",locus_without_gap[0].upper(),".",str(10000),"PASS",".","GT")
	elif len(locus_without_gap) == 2:
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("MT",str(i),".",locus_without_gap[0].upper(),locus_without_gap[1].upper(),str(10000),"PASS",".","GT")
	elif len(locus_without_gap) == 3:
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("MT",str(i),".",locus_without_gap[0].upper(),locus_without_gap[1].upper()+","+locus_without_gap[2].upper(),str(10000),"PASS",".","GT")
	elif len(locus_without_gap) == 4:
		line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format("MT",str(i),".",locus_without_gap[0].upper(),locus_without_gap[1].upper()+","+locus_without_gap[2].upper()+","+locus_without_gap[3].upper(),str(10000),"PASS",".","GT")
	
	for j in range(len(list(gene["NAME"]))):
		if gene["SEQUENCE"][j][i] == "-":
			line+="{}\t".format(".")
		else:
			line+="{}\t".format(str(locus_without_gap.index(gene["SEQUENCE"][j][i]))+"/.")
	
	with open("/home/labosea1/MitoZ/ANALYSIS/"+str(arg[1])+"/VCF/vcf_"+str(arg[1])+".vcf", "a") as file_object:
		file_object.write(line)
		file_object.write("\n")		

vcf=allel.read_vcf("/home/labosea1/MitoZ/ANALYSIS/"+str(arg[1])+"/VCF/vcf_"+str(arg[1])+".vcf",fields=['samples','calldata/GT'],log=sys.stdout)
samples = vcf['samples']

Li=[x for x in range(len(samples)) if 'Li' in samples[x]]
Mu=[x for x in range(len(samples)) if 'Mu' in samples[x]]
Fa=[x for x in range(len(samples)) if 'Fa' in samples[x]]
Ga=[x for x in range(len(samples)) if 'Ga' in samples[x]]
# Per poulations
subpops={	
	'all': range(len(vcf['samples'])),
	'Li': Li,
	'Mu': Mu,
	'Fa' : Fa,
	'Ga' : Ga,
	'Med' : Li+Mu,
	'Atl' : Fa+Ga,
	'Li-Fa' : Li+Fa,
	'Li-Ga' : Li+Ga,
	'Mu-Fa' : Mu+Fa,
	'Mu-Ga' : Mu+Ga,
	'Li-Mu-Fa' : Li+Mu+Fa,
	'Li-Mu-Ga' : Li+Mu+Ga,
	'Li-Fa-Ga' : Li+Fa+Ga,
	'Mu-Fa-Ga' : Mu+Fa+Ga
}

# All samples
gt = allel.GenotypeArray(vcf['calldata/GT']) #Convert vcf in GenotypeArray
h = gt.to_haplotypes #Convert in Haplotype format
ac=gt.count_alleles_subpops(subpops)


fst_hudson_Li_Mu, fst_hudson_se_Li_Mu, fst_hudson_vb_Li_Mu, _=allel.blockwise_hudson_fst(ac['Li'],ac['Mu'],blen=50)
print("Li-Mu Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Li_Mu,fst_hudson_se_Li_Mu))
print("-------------------------")
fst_hudson_Li_Fa, fst_hudson_se_Li_Fa, fst_hudson_vb_Li_Fa, _=allel.blockwise_hudson_fst(ac['Li'],ac['Fa'],blen=50)
print("Li-Fa Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Li_Fa,fst_hudson_se_Li_Fa))
print("-------------------------")
fst_hudson_Li_Ga, fst_hudson_se_Li_Ga, fst_hudson_vb_Li_Ga, _=allel.blockwise_hudson_fst(ac['Li'],ac['Ga'],blen=50)
print("Li-Ga Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Li_Ga,fst_hudson_se_Li_Ga))
print("-------------------------")
fst_hudson_Mu_Fa, fst_hudson_se_Mu_Fa, fst_hudson_vb_Mu_Fa, _=allel.blockwise_hudson_fst(ac['Mu'],ac['Fa'],blen=50)
print("Mu-Fa Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Mu_Fa,fst_hudson_se_Mu_Fa))
print("-------------------------")
fst_hudson_Mu_Ga, fst_hudson_se_Mu_Ga, fst_hudson_vb_Mu_Ga, _=allel.blockwise_hudson_fst(ac['Mu'],ac['Ga'],blen=50)
print("Mu-Ga Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Mu_Ga,fst_hudson_se_Mu_Ga))
print("-------------------------")
fst_hudson_Fa_Ga, fst_hudson_se_Fa_Ga, fst_hudson_vb_Fa_Ga, _=allel.blockwise_hudson_fst(ac['Fa'],ac['Ga'],blen=50)
print("Fa-Ga Hudson 's Fst: %.3f +/- %.3f" % (fst_hudson_Fa_Ga,fst_hudson_se_Fa_Ga))
print("-------------------------")

# Global Weir and Cockerham's Fst
fst_weirandcockerham_Li_Mu, fst_weirandcockerham_se_Li_Mu, fst_weirandcockerham_vb_Li_Mu, _=allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Li'],subpops['Mu']],blen=50,max_allele=1)
print("Li-Mu Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Li_Mu,fst_weirandcockerham_se_Li_Mu))
print("-------------------------")
fst_weirandcockerham_Li_Fa, fst_weirandcockerham_se_Li_Fa, fst_weirandcockerham_vb_Li_Fa, _= allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Li'],subpops['Fa']],blen=50,max_allele=1)
print("Li-Fa Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Li_Fa,fst_weirandcockerham_se_Li_Fa))
print("-------------------------")
fst_weirandcockerham_Li_Ga, fst_weirandcockerham_se_Li_Ga, fst_weirandcockerham_vb_Li_Ga, _= allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Li'],subpops['Ga']],blen=50,max_allele=1)
print("Li-Ga Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Li_Ga,fst_weirandcockerham_se_Li_Ga))
print("-------------------------")
fst_weirandcockerham_Mu_Fa, fst_weirandcockerham_se_Mu_Fa, fst_weirandcockerham_vb_Mu_Fa, _= allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Mu'],subpops['Fa']],blen=50,max_allele=1)
print("Mu-Fa Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Mu_Fa,fst_weirandcockerham_se_Mu_Fa))
print("-------------------------")
fst_weirandcockerham_Mu_Ga, fst_weirandcockerham_se_Mu_Ga, fst_weirandcockerham_vb_Mu_Ga, _= allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Mu'],subpops['Ga']],blen=50,max_allele=1)
print("Mu-Ga Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Mu_Ga,fst_weirandcockerham_se_Mu_Ga))
print("-------------------------")
fst_weirandcockerham_Fa_Ga, fst_weirandcockerham_se_Fa_Ga, fst_weirandcockerham_vb_Fa_Ga, _= allel.blockwise_weir_cockerham_fst(gt,subpops=[subpops['Fa'],subpops['Ga']],blen=50,max_allele=1)
print("Fa-Ga Weir & Cockerham 's Fst: %.3f +/- %.3f" % (fst_weirandcockerham_Fa_Ga,fst_weirandcockerham_se_Fa_Ga))
print("-------------------------")

tajimaD_all=allel.tajima_d(ac['all'])
print("All Tajima's D: %.3f" % (tajimaD_all))
tajimaD_li=allel.tajima_d(ac['Li'])
print("Li Tajima's D: %.3f" % (tajimaD_li))
tajimaD_mu=allel.tajima_d(ac['Mu'])
print("Li Tajima's D: %.3f" % (tajimaD_mu))
tajimaD_fa=allel.tajima_d(ac['Fa'])
print("Li Tajima's D: %.3f" % (tajimaD_fa))
tajimaD_ga=allel.tajima_d(ac['Ga'])
print("Li Tajima's D: %.3f" % (tajimaD_ga))

run_run = 0
if run_run == 1:
	is_bi_allelic=ac['all'].is_biallelic_01()[:]
	ac['Li'] = ac['Li'].compress(is_bi_allelic,axis=0)[:, :2]
	ac['Mu'] = ac['Mu'].compress(is_bi_allelic,axis=0)[:, :2]
	ac['Fa'] = ac['Fa'].compress(is_bi_allelic,axis=0)[:, :2]
	ac['Ga'] = ac['Ga'].compress(is_bi_allelic,axis=0)[:, :2]
	
	f3_Li_Mu_Fa, f3_se_Li_Mu_Fa, f3_z_Li_Mu_Fa,_,_=allel.average_patterson_f3(ac['Li'],ac['Mu'],ac['Fa'],100)
	print("Li;Mu,Fa f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Li_Mu_Fa,f3_se_Li_Mu_Fa,f3_z_Li_Mu_Fa))
	f3_Li_Mu_Ga, f3_se_Li_Mu_Ga, f3_z_Li_Mu_Ga,_,_=allel.average_patterson_f3(ac['Li'],ac['Mu'],ac['Ga'],100)
	print("Li;Mu,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Li_Mu_Ga,f3_se_Li_Mu_Ga,f3_z_Li_Mu_Ga))
	f3_Li_Fa_Ga, f3_se_Li_Fa_Ga, f3_z_Li_Fa_Ga,_,_=allel.average_patterson_f3(ac['Li'],ac['Fa'],ac['Ga'],100)
	print("Li;Fa,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Li_Fa_Ga,f3_se_Li_Fa_Ga,f3_z_Li_Fa_Ga))
	f3_Mu_Li_Fa, f3_se_Mu_Li_Fa, f3_z_Mu_Li_Fa,_,_=allel.average_patterson_f3(ac['Mu'],ac['Li'],ac['Fa'],100)
	print("Mu;Li,Fa f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Mu_Li_Fa,f3_se_Mu_Li_Fa,f3_z_Mu_Li_Fa))
	f3_Mu_Li_Ga, f3_se_Mu_Li_Ga, f3_z_Mu_Li_Ga,_,_=allel.average_patterson_f3(ac['Mu'],ac['Li'],ac['Ga'],100)
	print("Mu;Li,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Mu_Li_Ga,f3_se_Mu_Li_Ga,f3_z_Mu_Li_Ga))
	f3_Mu_Fa_Ga, f3_se_Mu_Fa_Ga, f3_z_Mu_Fa_Ga,_,_=allel.average_patterson_f3(ac['Mu'],ac['Fa'],ac['Ga'],100)
	print("Mu;Fa,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Mu_Fa_Ga,f3_se_Mu_Fa_Ga,f3_z_Mu_Fa_Ga))
	f3_Fa_Li_Mu, f3_se_Fa_Li_Mu, f3_z_Fa_Li_Mu,_,_=allel.average_patterson_f3(ac['Fa'],ac['Li'],ac['Mu'],100)
	print("Fa;Li,Mu f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Fa_Li_Mu,f3_se_Fa_Li_Mu,f3_z_Fa_Li_Mu))
	f3_Fa_Li_Ga, f3_se_Fa_Li_Ga, f3_z_Fa_Li_Ga,_,_=allel.average_patterson_f3(ac['Fa'],ac['Li'],ac['Ga'],100)
	print("Fa;Li,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Fa_Li_Ga,f3_se_Fa_Li_Ga,f3_z_Fa_Li_Ga))
	f3_Fa_Mu_Ga, f3_se_Fa_Mu_Ga, f3_z_Fa_Mu_Ga,_,_=allel.average_patterson_f3(ac['Fa'],ac['Mu'],ac['Ga'],100)
	print("Fa;Mu,Ga f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Fa_Mu_Ga,f3_se_Fa_Mu_Ga,f3_z_Fa_Mu_Ga))
	f3_Ga_Li_Mu, f3_se_Ga_Li_Mu, f3_z_Ga_Li_Mu,_,_=allel.average_patterson_f3(ac['Ga'],ac['Li'],ac['Mu'],100)
	print("Ga;Li,Mu f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Ga_Li_Mu,f3_se_Ga_Li_Mu,f3_z_Ga_Li_Mu))
	f3_Ga_Li_Fa, f3_se_Ga_Li_Fa, f3_z_Ga_Li_Fa,_,_=allel.average_patterson_f3(ac['Ga'],ac['Li'],ac['Fa'],100)
	print("Ga;Li,Fa f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Ga_Li_Fa,f3_se_Ga_Li_Fa,f3_z_Ga_Li_Fa))
	f3_Ga_Mu_Fa, f3_se_Ga_Mu_Fa, f3_z_Ga_Mu_Fa,_,_=allel.average_patterson_f3(ac['Ga'],ac['Mu'],ac['Fa'],100)
	print("Ga;Mu,Fa f3: %.3f +/- %.3f (Z = %.3f)" % (f3_Ga_Mu_Fa,f3_se_Ga_Mu_Fa,f3_z_Ga_Mu_Fa))


if str(arg[1])!="Afall":
	popgen_stats=pd.DataFrame(
	{
		"SP":arg[1],
		"GENE":arg[3],
		"LENGTH":length,
		"BIALLELIC":biallelic,
		"POLYALLELIC":polyallelic,
		"GAPS_ALL_INDIVIDUALS":gaps_all_individuals,
		"GAPS":gaps,
		"GC_CONTENT":gc_content,
		"PI_ALL":pi_all,
		"PI_LI":pi_li,
		"PI_MU":pi_mu,
		"PI_FA":pi_fa,
		"PI_GA":pi_ga,
		"DXY_LI_MU":dxy_li_mu,
		"DXY_LI_FA":dxy_li_fa,
		"DXY_LI_GA":dxy_li_ga,
		"DXY_MU_FA":dxy_mu_fa,
		"DXY_MU_GA":dxy_mu_ga,
		"DXY_FA_GA":dxy_fa_ga,		
		"DA_LI_MU":da_li_mu,
		"DA_LI_FA":da_li_fa,
		"DA_LI_GA":da_li_ga,
		"DA_MU_FA":da_mu_fa,
		"DA_MU_GA":da_mu_ga,
		"DA_FA_GA":da_fa_ga,
		"FST_HUDSON_LI_MU":fst_hudson_Li_Mu,
		"FST_HUDSON_LI_FA":fst_hudson_Li_Fa,
		"FST_HUDSON_LI_GA":fst_hudson_Li_Ga,
		"FST_HUDSON_MU_FA":fst_hudson_Mu_Fa,
		"FST_HUDSON_MU_GA":fst_hudson_Mu_Ga,
		"FST_HUDSON_FA_GA":fst_hudson_Fa_Ga,
		"FST_WEIRCOCKERHAM_LI_MU":fst_weirandcockerham_Li_Mu,
		"FST_WEIRCOCKERHAM_LI_FA":fst_weirandcockerham_Li_Fa,
		"FST_WEIRCOCKERHAM_LI_GA":fst_weirandcockerham_Li_Ga,
		"FST_WEIRCOCKERHAM_MU_FA":fst_weirandcockerham_Mu_Fa,
		"FST_WEIRCOCKERHAM_MU_GA":fst_weirandcockerham_Mu_Ga,
		"FST_WEIRCOCKERHAM_FA_GA":fst_weirandcockerham_Fa_Ga,
		"TAJIMAD_ALL":tajimaD_all,
		"TAJIMAD_LI":tajimaD_li,
		"TAJIMAD_MU":tajimaD_mu,
		"TAJIMAD_FA":tajimaD_fa,
		"TAJIMAD_GA":tajimaD_ga		
	},
	index=[0]
	)
else:
	popgen_stats=pd.DataFrame(
	{
		"SP":arg[1],
		"GENE":arg[3],
		"LENGTH":length,
		"BIALLELIC":biallelic,
		"POLYALLELIC":polyallelic,
		"GAPS_ALL_INDIVIDUALS":gaps_all_individuals,
		"GAPS":gaps,
		"GC_CONTENT":gc_content,
		"PI_ALL":pi_all,
		"PI_LI":pi_li,
		"PI_MU":"NA",
		"PI_FA":pi_fa,
		"PI_GA":pi_ga,
		"DXY_LI_MU":"NA",
		"DXY_LI_FA":dxy_li_fa,
		"DXY_LI_GA":dxy_li_ga,
		"DXY_MU_FA":"NA",
		"DXY_MU_GA":"NA",
		"DXY_FA_GA":dxy_fa_ga,		
		"DA_LI_MU":"NA",
		"DA_LI_FA":da_li_fa,
		"DA_LI_GA":da_li_ga,
		"DA_MU_FA":"NA",
		"DA_MU_GA":"NA",
		"DA_FA_GA":da_fa_ga,
		"FST_HUDSON_LI_MU":fst_hudson_Li_Mu,
		"FST_HUDSON_LI_FA":fst_hudson_Li_Fa,
		"FST_HUDSON_LI_GA":fst_hudson_Li_Ga,
		"FST_HUDSON_MU_FA":fst_hudson_Mu_Fa,
		"FST_HUDSON_MU_GA":fst_hudson_Mu_Ga,
		"FST_HUDSON_FA_GA":fst_hudson_Fa_Ga,
		"FST_WEIRCOCKERHAM_LI_MU":fst_weirandcockerham_Li_Mu,
		"FST_WEIRCOCKERHAM_LI_FA":fst_weirandcockerham_Li_Fa,
		"FST_WEIRCOCKERHAM_LI_GA":fst_weirandcockerham_Li_Ga,
		"FST_WEIRCOCKERHAM_MU_FA":fst_weirandcockerham_Mu_Fa,
		"FST_WEIRCOCKERHAM_MU_GA":fst_weirandcockerham_Mu_Ga,
		"FST_WEIRCOCKERHAM_FA_GA":fst_weirandcockerham_Fa_Ga,
		"TAJIMAD_ALL":tajimaD_all,
		"TAJIMAD_LI":tajimaD_li,
		"TAJIMAD_MU":tajimaD_mu,
		"TAJIMAD_FA":tajimaD_fa,
		"TAJIMAD_GA":tajimaD_ga		
	},
	index=[0]
	)	
popgen_stats.to_csv("/home/labosea1/MitoZ/ANALYSIS/"+str(arg[1])+"/Pop_Gen/popgen_stats_"+str(arg[3])+".csv",index=False)



















