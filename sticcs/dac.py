#!/usr/bin/env python

import argparse, sys, gzip

import cyvcf2
import numpy as np

def get_allele_counts(geno, n_alleles):
    missing = np.array([-1]*n_alleles)
    try: return np.bincount(geno, minlength = n_alleles)
    except: return missing

def is_informative(dercounts, minmax):
    dersum = dercounts[dercounts>0].sum()
    if minmax[0] and minmax[0] > dersum: return False
    if minmax[1] and minmax[1] < dersum: return False
    return True

def get_anc_from_outgroup(outgroup_allele_counts):
    usable = np.where(outgroup_allele_counts.min(axis = 1) != -1)
    if len(usable[0]) >= 1:
        outgroup_allele_totals = outgroup_allele_counts[usable].sum(axis=0)
        anc = np.where(outgroup_allele_totals > 0)[0]
        if len(anc) == 1: return anc[0]
    return None

def get_anc_from_field(variant):
    ancestral = variant.INFO.get("AA", ".")  # "." means unknown
    # some VCFs (e.g. from 1000G) have many values in the AA field: take the 1st
    ancestral = ancestral.split("|")[0].upper()
    if ancestral == "." or ancestral == "": return None
    return ancestral


def allele_counts_to_dacs(allele_counts, anc, minmax):
    n_alleles = allele_counts.shape[1]
    
    assert n_alleles > 1 and anc < n_alleles
    
    ders = [i for i in range(n_alleles) if i != anc]
    
    site_dercounts = allele_counts[:,ders]
    
    if n_alleles > 2:
        #set individuals with more than one derived allele to missing (These are unusable by sticcs)
        site_dercounts[(site_dercounts > 0).sum(axis=1) > 1] = -1
        
        #for each derived allele, we cannot consider individuals that carry any other derived alleles, so set these to missing
        for i in range(n_alleles-1):
            der_present = site_dercounts[:,i] > 0
            if der_present.sum() > 0:
                for other in [j for j in range(n_alleles-1) if j!=i]:
                    site_dercounts[der_present, other] = -1
    
    if minmax:
        site_dercounts = site_dercounts[:,np.apply_along_axis(is_informative, 0, site_dercounts, minmax)]
    
    return site_dercounts.transpose()


#extract derived counts from a haploid matrix. This is useful for processing simulated data
#assumes 0 is ancestral
def get_dac_from_haploid_matrix(mat, positions, ploidies):
    
    counts_list = []
    new_positions = []
    
    partitions = np.cumsum(ploidies)[:-1]
    
    total=ploidies.sum()
    
    n = len(positions)
    assert n == mat.shape[0]
    
    for i in range(n):
        n_alleles = mat[i].max() + 1
        if n_alleles > 1:
            genos = np.split(mat[i], partitions)
            allele_counts = np.array([get_allele_counts(geno, n_alleles) for geno in genos])
            site_dercounts = allele_counts_to_dacs(allele_counts, anc=0, minmax=(2, total-1,))
            counts_list.append(site_dercounts)
            #for each derived allele that we have a counts list for, we add the position to the list
            for j in range(site_dercounts.shape[0]):
                new_positions.append(positions[i])
    
    return (np.vstack(counts_list), np.array(new_positions))


def make_vcf_with_DC_field(vcf, out_name, outgroups=None, use_REF_as_anc=False):
    
    if outgroups is None and not use_REF_as_anc:
        try: vcf.get_header_type("AA")
        except: raise Exception("The VCF has no ancestral allele (AA) field. To infer this you must specify outgroups.")
    
    outgroup_indices = [vcf.samples.index(sample) for sample in outgroups] if outgroups else None

    #add necessary header line
    vcf.add_format_to_header({"ID":"DC", "Number":1, "Type":"Integer", "Description":"Derived allele count"})
    
    out = cyvcf2.Writer(out_name, vcf)
    
    out.write_header()
    
    for variant in vcf:
        n_alleles = len(variant.ALT) + 1
        
        if n_alleles < 2: continue
        
        allele_counts = np.array([get_allele_counts(geno[:2], n_alleles) for geno in variant.genotypes])
        
        if use_REF_as_anc:
            anc=0
        elif outgroup_indices:
            outgroup_genos = [gt[:-1] for i,gt in enumerate(variant.genotypes) if i in outgroup_indices]
            outgroup_allele_counts = np.array([get_allele_counts(geno, n_alleles) for geno in outgroup_genos])
            anc = get_anc_from_outgroup(outgroup_allele_counts)
        else:
            anc = get_anc_from_field(variant)
        
        if anc is None:
            continue
        
        for dac in allele_counts_to_dacs(allele_counts, anc, minmax=(2, None,)):
            variant.set_format("DC", dac)
            out.write_record(variant)
        
    out.close()

#functions for reading files that already have derived count info
#yields derived couns for one chromosome at a time
def parse_vcf_with_DC_field(vcf):
    
    positions = []
    dacs = []
    last_chrom = None
    
    for variant in vcf:
        if last_chrom is None: last_chrom = variant.CHROM
        
        dac = variant.format("DC").transpose()[0]
        
        if variant.CHROM == last_chrom:
            positions.append(variant.POS)
            dacs.append(dac)
        
        else:
            yield (last_chrom, np.array(positions, dtype=int), np.array(dacs, dtype=int))
            last_chrom = variant.CHROM
            positions = [variant.POS]
            dacs = [dac]
    
    yield (last_chrom, np.array(positions, dtype=int), np.array(dacs, dtype=int))


###############################################################################################################
if __name__ == "__main__":
    
    ### parse arguments

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in_file", help="Input vcf file", action = "store", default="-")
    parser.add_argument("-o", "--out_file", help="Output vcf file", action = "store", default="-")
    parser.add_argument("--outgroups", help="sample name(s) for outgroup(s)", action='store', nargs="+")
    parser.add_argument("--use_REF_as_anc", help="Use REF allele as ancestral.", action = "store_true")
    
    args = parser.parse_args()
    
    vcf = cyvcf2.VCF(args.in_file)
    
    make_vcf_with_DC_field(vcf, args.out_file, outgroups=args.outgroups, use_REF_as_anc=args.use_REF_as_anc)
