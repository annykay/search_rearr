import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy import stats
from scipy.stats import nbinom
from patsy import dmatrices
import re 
import random
import Bio
from Bio import SeqIO
import sys

def read_pairs(genome):
    genome = genome.replace('-', '_')
    first_filepath = 'SourceData/results_' + genome + \
                     '/gmapped_parsed_sorted_chunks/coverages/chimersFirst_inters.cov'
    second_filepath = 'SourceData/results_' + genome + \
                      '/gmapped_parsed_sorted_chunks/coverages/chimersSecond_inters.cov'
    print(first_filepath)
    first_mate = pd.read_csv(first_filepath, sep = '\t', index_col =False, 
                       header=None)
    second_mate = pd.read_csv(second_filepath, sep = '\t', index_col =False, 
                       header=None)
    first_mate = first_mate[[0,1,3]]
    second_mate = second_mate[[0,1,3]]

    first_mate.columns = ['chrom1', '3end1', 'restr1']
    second_mate.columns = ['chrom2', '3end2', 'restr2']
    pairs = pd.merge(first_mate, second_mate, left_index=True, right_index=True)
    pairs = pairs[pairs['restr1'] == 0]
    pairs = pairs[pairs['restr2'] == 0]
    return pairs

def read_chrom_size(genome):
    chrom_size_file = 'SourceData/indexes/' + genome + '_scaffolds.chrom.sizes'
    chrom_size = pd.read_csv(chrom_size_file, sep = '\t', header=None)
    chrom_size.columns = ['chrom', 'len']
    chrom_size['cum_len'] = np.cumsum(chrom_size['len']) - chrom_size['len'][0]
    return chrom_size

def filter_scaff(pairs, scaff):
    tmp_float = [i == scaff or j == scaff for i, j in zip(pairs['chrom1'], pairs['chrom2'])]
    tmp_idx = [i for i, x in enumerate(tmp_float) if x]
    return pairs.iloc[tmp_idx]

def filter_main_diag(pairs_one_scaff, chrom_size, width = 10000):
    pairs_one_scaff = pd.merge(pairs_one_scaff, chrom_size, left_on = 'chrom1', right_on = 'chrom')
    pairs_one_scaff['cum_3end1'] = pairs_one_scaff['3end1'] + pairs_one_scaff['cum_len']
    pairs_one_scaff = pairs_one_scaff.drop(['chrom', 'len', 'cum_len'], axis = 1)

    pairs_one_scaff = pd.merge(pairs_one_scaff, chrom_size, left_on = 'chrom2', right_on = 'chrom')
    pairs_one_scaff['cum_3end2'] = pairs_one_scaff['3end2'] + pairs_one_scaff['cum_len']
    pairs_one_scaff = pairs_one_scaff.drop(['chrom', 'len', 'cum_len'], axis = 1)

    return pairs_one_scaff[abs(pairs_one_scaff['cum_3end2'] - pairs_one_scaff['cum_3end1']) > 10000]

def find_imp_restr_regions(pairs_one_scaff, chrom_size, scaff, bin_size = 10):
    scaff_len = int(chrom_size[chrom_size['chrom'] == scaff]['len'])

    for_test = np.concatenate([pairs_one_scaff[pairs_one_scaff['chrom1'] == scaff]['3end1'],
                               pairs_one_scaff[pairs_one_scaff['chrom2'] == scaff]['3end2']])
    hist = np.histogram(for_test, 
                     bins = scaff_len // bin_size)
    hist_z = (hist[0] - np.mean(hist[0])) / np.std(hist[0])
    tre = np.percentile(hist_z, 99.99)
    tmp_idx = [i for i, x in enumerate(hist_z > tre) if x]
    return hist[1][tmp_idx]

def find_sequences(genome_file, regions, bin_size, scaff):
    with open(genome_file) as handle:
        sequences = []
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == scaff:
                for i in regions:
                    sequences.append(record.seq[int(i):(int(i) + 10)])
    return sequences

def check_sequences(sequences):
    result = 0
    for j in sequences:
        if re.findall(r"[A,C,T]ATC", str(j)):
            result += 1
        elif re.findall(r"G[C,T,G]TC", str(j)):
            result += 1
        elif re.findall(r"GA[A,C,G]C", str(j)):
            result += 1
        elif re.findall(r"GAT[A,T,G]", str(j)):
            result += 1
    return result

def make_implicit_sites(sequences, starts, scaff):
    mutation_points = []
    changes = []
    references = []
    restr_sites = []
    for i, seq in enumerate(sequences):
        if re.findall(r"[A,C,T]ATC", str(seq)):
            span = re.search(r"[A,C,T]ATC", str(seq)).span()
            mut_pos = starts[i] + span[0]
            ref = seq[span[0]]
            change = 'G'
            restr = starts[i] + span[0]

        elif re.findall(r"G[C,T,G]TC", str(seq)):
            span = re.search(r"G[C,T,G]TC", str(seq)).span()
            mut_pos = starts[i] + span[0] + 1
            ref = seq[span[0] + 1]
            change = 'A'
            restr = starts[i] + span[0]

        elif re.findall(r"GA[A,C,G]C", str(seq)):
            span = re.search(r"GA[A,C,G]C", str(seq)).span()
            mut_pos = starts[i] + span[0] + 2
            ref = seq[span[0] + 2]
            change = 'T'
            restr = starts[i] + span[0]

        elif re.findall(r"GAT[A,T,G]", str(seq)):
            span = re.search(r"GAT[A,T,G]", str(seq)).span()
            mut_pos = starts[i] + span[0] + 3
            ref = seq[span[0] + 3]
            change = 'C'
            restr = starts[i] + span[0]
        else: 
            mut_pos = -1
            ref = '.'
            change = '.'
            restr = -1
        mutation_points.append(mut_pos)
        changes.append(change)
        references.append(ref)
        restr_sites.append(restr)

    mutations = {'CHROM': scaff,
                 'POS': mutation_points,
                 'REF': references,
                 'ALT': changes}
    mutations = pd.DataFrame(data = mutations)
    mutations = mutations[mutations['REF'] != '.']
    
    restr_sites = {'CHROM': scaff,
                   'START': restr_sites}
    restr_sites = pd.DataFrame(data = restr_sites)
    restr_sites = restr_sites[restr_sites['START'] != -1]
    restr_sites['END'] = restr_sites['START'] + 4
    return restr_sites, mutations
def find_significance(genome_file, bin_size, scaff, scaff_len, n_found, n_confirm, n_tests = 20):
    perc = []
    for i in range(n_tests):
        rand_ind = random.sample(range(scaff_len), k = n_found)
        sequences = find_sequences(genome_file, rand_ind, bin_size, scaff)

        result = check_sequences(sequences)
        perc.append(result/n_found)
        
    return sum([i > n_confirm/n_found for i in perc])/n_tests * 100
    
    
    
def find_implicit_restriction_sites(pairs, scaff, chrom_size, genome_file, main_diag_width = 10000,
                                   bin_size = 10):
    scaff_len = int(chrom_size[chrom_size['chrom'] == scaff]['len'])
    #first, filter only interested chromosome
    pairs_one_scaff = filter_scaff(pairs, scaff)
    
    #Than, filter out elements on main diagonal
    pairs_one_scaff = filter_main_diag(pairs_one_scaff, chrom_size, main_diag_width)
    
    #find bins, where an outstanding amount of 3'wnd chimeric reads goes 
    impl_restr_region_starts = find_imp_restr_regions(pairs_one_scaff, chrom_size, 
                                                      scaff, bin_size)
    
    #calculate, how many of them have implicit restriction sites
    sequences = find_sequences(genome_file, impl_restr_region_starts, bin_size, scaff)
    n_implicit = check_sequences(sequences)
    
    #make a dataframe with restriction sites starts and point mutations lead to it 
    implicit_restriction_sites, mutations = make_implicit_sites(sequences, impl_restr_region_starts, scaff)

    #check significance of the result
    signif = find_significance(genome_file, bin_size, scaff, scaff_len, len(sequences), n_implicit)
    print(f'Significance: {signif} %')
    return mutations, implicit_restriction_sites
 
def rearrangments_search(pairs, chrom_size, chrom_1, chrom_2, bin_size, sign):
    chrom_length_1 = int(chrom_size[chrom_size['chrom'] == chrom_1]['len'])
    chrom_length_2 = int(chrom_size[chrom_size['chrom'] == chrom_2]['len'])


    n_bins = int(max(chrom_length_1, chrom_length_2) / bin_size)
    pairs_part = pairs[pairs['chrom1'] == chrom_1]
    pairs_part = pairs_part[pairs_part['chrom2'] == chrom_2]
    freqs, first, second = np.histogram2d(pairs_part['3end1'], pairs_part['3end2'], 
                                          bins = n_bins)
    freqs_1d = freqs.reshape(-1)
    tre = np.percentile(freqs_1d, sign)
    first_starts = []
    second_starts = []
    for i, first_start in enumerate(first):
        for j, second_start in enumerate(second):
            if i < len(freqs) and j < len(freqs[0]) and freqs[i][j] > tre:
                first_starts.append(int(first_start))
                second_starts.append(int(second_start))

    rearrangments = {'start1': first_starts, 'start2': second_starts}
    rearrangments = pd.DataFrame(rearrangments)

    rearrangments['chrom1'] = chrom_1
    rearrangments['chrom2'] = chrom_2
    rearrangments['end1'] = rearrangments['start1'] + bin_size
    rearrangments['end2'] = rearrangments['start2'] + bin_size
    rearrangments = rearrangments[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
    return rearrangments

def main():
    genome = sys.argv[1]
    genome_file = 'SourceData/indexes/' + genome + '_scaffolds.fna'
    chrom_size = read_chrom_size(genome)
    pairs = read_pairs(genome)
    chrom_names = chrom_size['chrom'][0:11]
    bin_size = 1000
    for scaff in chrom_names:
        mutations, imlicit_restr_sites = find_implicit_restriction_sites(pairs, scaff, chrom_size, genome_file)
        filename = 'Results/mutations/s' + genome + '_' + str(bin_size) + '_' + scaff + '.csv'
        mutations.to_csv(filename, sep = '\t', index = False)
        filename = 'Results/imlicit_restr_sites/s' + genome + '_' + str(bin_size) + '_' + scaff + '.csv'
        imlicit_restr_sites.to_csv(filename, sep = '\t', index = False)
    for i, chrom_1 in enumerate(chrom_names):
        for chrom_2 in chrom_names[i:]:
            df = rearrangments_search(pairs, chrom_size, chrom_1, chrom_2, bin_size, 99.9999)
            filename = 'Results/my_rearrangments/' + genome + '_' + str(bin_size) + '_' + chrom_1 + '_' + chrom_2 + '.csv'
            df.to_csv(filename, sep = '\t', index = False)
    
    

if __name__ == "__main__":
    main()
