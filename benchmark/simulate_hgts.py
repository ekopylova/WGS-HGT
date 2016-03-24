# ----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

#
# This script allows to simulate horizontal gene transfer by combining genes
# from raw nucleotide genomes.
#
# Input: the inputs are two genomes, donor and recipient (in GenBank format or
#        in single line format for artificial genomes by Azad & Lawrence) and
#        options regarding the number of HGTs to simulate and the types of
#        HGTs allowed (orthologous replacement and novel gene acquisition
#        supported)
#
# Algorithm: all protein coding sequences (cds) for both genomes are extracted
#            using scikit-bio. For orthologous gene replacement, orthologous
#            genes are determined using OrthoFinder. A python function
#            combines genes from the two genomes (orthologous and novel).
#
# Output:    the resulting genomes in raw nucleotide and protein coding
#            sequences FASTA formats and a log file showing the nucleotide
#            positions of spiked genes
#
# Dependencies: this script requires the tools FeatureExtract and OrthoFinder
#
# Additional information:
#   1. Emms, D. and Kelly, S. (2015). OrthoFinder: solving fundamental biases
#      in whole genome comparisons dramatically improves orthogroup inference
#      accuracy, GenomeBiology, 16:157
#   2. Azad, R. and Lawrence, J.G. (2005). Use of Artificial Genomes in
#      Assessing Methods for Atypical Gene Detection, PLoS Computational
#      Biology, doi: 10.1371/journal.pcbi.0010056
#

import sys
import click
from os import makedirs
from os.path import exists, basename, join, splitext
import subprocess
import time
import glob
import random
from operator import itemgetter

from skbio import Sequence, RNA, GeneticCode


def extract_azad_lawrence(artificial_genome_fp,
                          artificial_annotation_fp,
                          genome_id='seq'):
    """Extract protein coding sequences from artificial genomes.

    Parameters
    ----------
    artificial_genome_fp: string
        file path to artificial genome
    artificial_annotation_fp: string
        file path to annotation for artificial genome
    genome_id: string, optional
        unique genome ID

    Returns
    -------
    seq: skbio.sequence.Sequence
        Sequence object for raw nucleotide genome
    genes: dictionary
        a dictionary of genes (CDS) and their info, with the key being the
        protein IDs and the value being a 4-element list including the
        translated sequence, the start and end positions in the genome

    Notes
    -----
    Artificial genome must be one of those provided by (Azad & Lawrence, 2005)
    """
    genes = {}
    with open(artificial_genome_fp, 'U') as artificial_genome_f:
        seq = Sequence(artificial_genome_f.readline())
    with open(artificial_annotation_fp, 'U') as artificial_annotation_f:
        gene_positions = [line.strip().split()
                          for line in artificial_annotation_f]
    for pos in gene_positions:
        strand = pos[0]
        gene_start = int(pos[1])
        gene_end = int(pos[2])
        protein_id = "%s_%s_%s" % (genome_id, gene_start, gene_end)
        if protein_id not in genes:
            gene_nucl = str(seq[gene_start-1:gene_end-1])
            gene_nucl_rna = RNA(gene_nucl.replace('T', 'U'))
            gc = GeneticCode.from_ncbi(11)
            reading_frames = [1,2,3,-1,-2,-3]
            gene_trans_list = []
            for frame in reading_frames:
                # only record trimmed translations for the CDS (must contain
                # a start and stop codon)
                try:
                    gene_trans = gc.translate(
                        gene_nucl_rna,
                        reading_frame=frame,
                        start='require',
                        stop='require')
                except ValueError:
                    gene_trans = None
                if gene_trans:
                    gene_trans_list.append(gene_trans)
            # choose the longest translated sequence with start and stop
            # codons
            if gene_trans_list:
                gene_trans = max(gene_trans_list, key=len)
                genes[protein_id] = [str(gene_trans),
                                     gene_start-1,
                                     gene_end-1,
                                     strand]
        else:
            raise KeyError("%s already exists in dictionary" % protein_id)
    return seq, genes


def extract_genbank(genbank_fp, verbose=False):
    """Extract protein coding sequences from GenBank record.
    Parameters
    ----------
    genbank_fp: string
        file path to genome in GenBank format

    Returns
    -------
    seq: skbio.sequence.Sequence
        Sequence object
    genes: dictionary
        a dictionary of genes (CDS) and their info, with the key being the
        protein IDs and the value being a 4-element list including the
        translated sequence, the start and end positions in the genome
    """
    genes = {}
    if verbose:
        sys.stdout.write("\tParse GenBank record ...\n")
    seq = Sequence.read(genbank_fp, format='genbank')
    if verbose:
        sys.stdout.write("\t\tDone.\n")
    for feature in seq.metadata['FEATURES']:
        if feature['type_'] == 'CDS':
            protein_id = feature['protein_id']
            translation = feature['translation']
            strand = '+'
            if feature['rc_']:
                strand = '-'
            col = feature['index_']
            loc = seq.positional_metadata[col]
            start_pos = loc[loc == True].index[0]
            end_pos = loc[loc == True].index[-1]
            if protein_id not in genes:
                genes[protein_id.replace("\"", "")] = [
                    translation.replace(" ", "").replace("\"", ""),
                    start_pos, end_pos, strand]
            else:
                raise KeyError("%s already exists in dictionary" % protein_id)
    return seq, genes 


def launch_orthofinder(proteomes_dir, threads):
    """Launch OrthoFinder to report orthologous gene groups for two genomes.

    Parameters
    ----------
    proteomes_dir: string
        directory path storing FASTA coding sequences for complete genomes
    threads: integer
        number of threads to use
    """
    orthofinder_command = ["orthofinder.py",
                           "-f", proteomes_dir,
                           "-t", str(threads)]
    proc = subprocess.Popen(orthofinder_command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr


def parse_orthofinder(results_dir):
    """Parse the output files of OrthoFinder for orthologous genes.

    Parameters
    ----------
    results_dir: string
        OrthoFinder results directory

    Returns
    -------
    species_ids: dictionary
    sequence_ids: dictionary
    orthologous_groups: list of lists
    """
    species_ids = {}
    sequence_ids = {}
    orthologous_groups = []
    with open(
        glob.glob(
            join(
                results_dir,
                "clusters_OrthoFinder_*_id_pairs.txt"))[0], 'U') as id_pairs_f:
        # skip header lines
        for _ in xrange(7):
            next(id_pairs_f)
        # parse orthologous family groups
        for line in id_pairs_f:
            line = line.strip().split()
            # include only families with at least 2 orthologs
            if len(line[1:-1]) > 1:
                orthologous_groups.append(line[1:-1])
    with open(join(results_dir, "SpeciesIDs.txt"), 'U') as species_ids_f:
        for line in species_ids_f:
            line = line.strip().split()
            species_ids[line[0].split(':')[0]] = line[1]
    with open(join(results_dir, "SequenceIDs.txt"), 'U') as sequence_ids_f:
        for line in sequence_ids_f:
            line = line.strip().split()
            sequence_ids[line[0].split(':')[0]] = line[1]

    return species_ids, sequence_ids, orthologous_groups


def simulate_orthologous_rep(genes_donor,
                             seq_donor,
                             genes_recip,
                             seq_recip,
                             sequence_ids,
                             orthologous_groups,
                             orthologous_rep_prob,
                             percentage_hgts,
                             log_f):
    """Simulate orthologous replacement HGT.

    Notes
    -----
    Algorithm: Using list of orthologous genes between donor and recipient
    genomes 

    Parameters
    ----------
    genes_donor: dictionary
    seq_donor: skbio.sequence.Sequence
    genes_recip: dictionary
    seq_recip: skbio.sequence.Sequence
    sequence_ids: dictionary
    orthologous_groups: list of lists
    orthologous_rep_prob: float
    percentage_hgts: float
    log_f: file descriptor

    Returns
    -------
    seq_recip: skbio.sequence.Sequence
        recipient genome sequence with HGTs
    """
    # number of HGTs to simulate
    num_hgts = int(percentage_hgts*orthologous_rep_prob*len(genes_recip))
    if num_hgts < 1:
        num_hgts = 1
    num_orthogroups = len(orthologous_groups)
    # 'num_hgts' random number of indexes for orthologous_groups list
    idx = random.sample(xrange(0, num_orthogroups), num_hgts)
    log_f.write("#type\tdonor\tstart\tend\trecipient\tnew label "
                "recipient\tstart\tend\tstrand\n")
    for x in xrange(0, num_hgts):
        orthogroup = orthologous_groups[idx[x]]
        substitute_genes = ['*', '*']
        # randomly select two orthologous genes from the same family
        # representing the donor and recipient genomes 
        while '*' in substitute_genes:
            idx2 = random.randrange(0, len(orthogroup))
            gene = orthogroup[idx2]
            if (gene.startswith('0') and substitute_genes[0] == '*'):
                substitute_genes[0] = sequence_ids[gene]
            elif (gene.startswith('1') and substitute_genes[1] == '*'):
                substitute_genes[1] = sequence_ids[gene]
        # match donor and recipient gene labels to results output by
        # OrthoFinder (in sequence_ids)
        gene_donor_label = None
        gene_recip_label = None
        if substitute_genes[0] in genes_donor:
            gene_donor_label = substitute_genes[0]
            gene_recip_label = substitute_genes[1]
        elif substitute_genes[1] in genes_donor:
            gene_donor_label = substitute_genes[1]
            gene_recip_label = substitute_genes[0]
        else:
            raise ValueError("Gene %s and %s are not in donor genome" %
                (substitute_genes[0], substitute_genes[1]))
        # rename recipient orthologous gene to donor's
        hgt_gene = "%s_hgt_o" % gene_donor_label
        genes_recip[hgt_gene] = genes_recip.pop(gene_recip_label)
        # replace recipient gene (translated sequence) with donor's
        genes_recip[hgt_gene][0] = genes_donor[gene_donor_label][0]
        # update end position of HGT gene (as it can be shorter/longer than
        # the recipient gene replaced)
        genes_recip[hgt_gene][2] =\
            genes_recip[hgt_gene][1] + len(genes_recip[hgt_gene][0])*3
        # replace recipient gene (nucleotide format) with donor's
        start_pos_recip, end_pos_recip, strand_recip =\
            genes_recip[hgt_gene][1:]
        start_pos_donor, end_pos_donor, strand_donor =\
            genes_donor[gene_donor_label][1:]
        seq_recip = Sequence(str(seq_recip[:start_pos_recip]) +\
            str(seq_donor[start_pos_donor:end_pos_donor]) +\
            str(seq_recip[end_pos_recip:]))
        if strand_recip != strand_donor:
            genes_recip[hgt_gene][3] =\
            genes_donor[gene_donor_label][3]
        # write HGTs to log file
        log_f.write("o\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            gene_donor_label,
            start_pos_donor,
            end_pos_donor,
            gene_recip_label,
            hgt_gene,
            start_pos_recip,
            end_pos_recip,
            strand_donor))
    return seq_recip


def simulate_novel_acq(genes_donor,
                       seq_donor,
                       genes_recip,
                       seq_recip,
                       orthologous_rep_prob,
                       percentage_hgts,
                       log_f):
    """Simulate novel gene acquisition HGT.

    Parameters
    ----------
    genes_donor: dictionary
    seq_donor: skbio.sequence.Sequence
    genes_recip: dictionary
    seq_recip: skbio.sequence.Sequence
    orthologous_rep_prob: float
    percentage_hgts: float 
    log_f: file descriptor

    Returns
    -------
    seq_recip: skbio.sequence.Sequence
        recipient genome sequence with HGTs

    Notes
    -----
    Algorithm:
        1. choose random location in recipient genome where to insert a gene
           (randomly chosen from list of donor genes)
        2. use gene (recipient) positioning array to locate an open region
           (that doesn't include an existing gene) near the random location
           to insert the new gene (we want to avoid gene overlap so that
           compositional methods can clearly pick out individual coding
           genes)
        3. insert new gene, record existance in gene positioning array
    """
    num_hgts = int(percentage_hgts*(1-orthologous_rep_prob)*len(genes_recip))
    if num_hgts < 1:
        num_hgts = 1
    # create recipient genome gene positioning array
    gene_positions =\
        [(genes_recip[gene][1], genes_recip[gene][2]) for gene in genes_recip]
    # add start and end positions of recipient genome to allow for HGTs
    # simulated before the first and after the last existing gene
    gene_positions.append((0,0))
    gene_positions.append((len(seq_recip), len(seq_recip)))
    # sort array for gene positions in ascending order
    gene_positions_s = sorted(gene_positions, key=itemgetter(0))
    random_used = []
    end_simulation = False
    # select a random list of positions where to insert the new gene
    idx = random.sample(xrange(0, len(gene_positions_s)-1), num_hgts)
    gene_donor_labels = random.sample(genes_donor.keys(), num_hgts)
    log_f.write("#type\tdonor\tstart\tend\trecipient\tstart\t"
                "end\tstrand\n")
    # begin simulation
    for x in xrange(0, num_hgts):
        # select random donor gene (for HGT)
        gene_donor_label = gene_donor_labels[x]
        idx_recip = gene_positions_s[idx[x]][1] + 1
        # beginning from valid position for inserting new gene, check whether
        # the length can fit without overlapping with existing gene, otherwise
        # search for next valid position
        for y in xrange(idx[x], len(gene_positions_s)-1):
            if idx_recip + len(genes_donor[gene_donor_label][0])*3 <\
                    gene_positions_s[y+1][0]:
                idx_end = idx_recip + len(genes_donor[gene_donor_label][0])*3
                # insert gene (protein)
                hgt_gene = "%s_hgt_n" % gene_donor_label
                genes_recip[hgt_gene] =\
                    [genes_donor[gene_donor_label][0], idx_recip, idx_end,
                     genes_donor[gene_donor_label][3]]
                # insert gene (nucleotide)
                seq_recip = Sequence(str(seq_recip[:idx_recip]) +\
                    str(seq_donor[genes_donor[gene_donor_label][1]:\
                        genes_donor[gene_donor_label][2]]) +\
                    str(seq_recip[idx_recip:]))
                # write HGTs to log file
                log_f.write(
                    "n\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %\
                    (gene_donor_label,
                     genes_donor[gene_donor_label][1],
                     genes_donor[gene_donor_label][2],
                     hgt_gene,
                     idx_recip,
                     idx_end,
                     genes_donor[gene_donor_label][3]))
                break
            # try next open region
            idx_recip = gene_positions_s[y+1][1] + 1

    return seq_recip


def write_results(genes_donor,
                  donor_genome_fp,
                  genes_recip,
                  recip_genome_fp,
                  seq_donor,
                  seq_recip,
                  output_dir):
    """Write donor and recipient genomes to FASTA files (nucl and protein).

    Parameters
    ----------
    genes_donor: dictionary
    seq_donor: skbio.sequence.Sequence
        Sequence object for donor genome
    genes_recip: dictionary
    seq_recip: skbio.sequence.Sequence
        Sequence object for recipient genome
    """
    # output dir for simulated results
    simulated_dir = join(output_dir, "simulated")
    if not exists(simulated_dir):
        makedirs(simulated_dir)
    donor_genome_aa_fp = join(
        simulated_dir, "%s.faa" % basename(splitext(donor_genome_fp)[0]))
    recip_genome_aa_fp = join(
        simulated_dir, "%s.faa" % basename(splitext(recip_genome_fp)[0]))
    with open(donor_genome_aa_fp, 'w') as donor_genome_aa_f:
        for gene in genes_donor:
            donor_genome_aa_f.write(
                ">%s\n%s\n" % (gene, genes_donor[gene][0]))
    with open(recip_genome_aa_fp, 'w') as recip_genome_aa_f:
        for gene in genes_recip:
            recip_genome_aa_f.write(
                ">%s\n%s\n" % (gene, genes_recip[gene][0]))
    donor_genome_nucl_fp = join(
        simulated_dir, "%s.fna" % basename(splitext(donor_genome_fp)[0]))
    recip_genome_nucl_fp = join(
        simulated_dir, "%s.fna" % basename(splitext(recip_genome_fp)[0]))
    seq_donor.write(donor_genome_nucl_fp, format='fasta')
    seq_recip.write(recip_genome_nucl_fp, format='fasta')

    return (donor_genome_nucl_fp, donor_genome_aa_fp, recip_genome_nucl_fp,
        recip_genome_aa_fp)

def simulate_hgts(seq_donor,
                  genes_donor,
                  seq_recip,
                  genes_recip,
                  donor_genome_fp,
                  recip_genome_fp,
                  output_dir,
                  percentage_hgts,
                  orthologous_rep_prob,
                  log_f,
                  threads=1,
                  verbose=False):
    """Simulate orthologous replacement and novel gene acquisition HGTs.

    Parameters
    ----------
    seq_donor: skbio.sequence.Sequence
        Sequence object for donor genome
    genes_donor: dictionary
        a dictionary of genes (CDS) and their info, with the key being the
        protein IDs and the value being a 5-element list including the
        translated sequence, the raw nucleotide sequence, the start and end
        positions in the donor genome
    seq_recip: skbio.sequence.Sequence
        Sequence object for recipient genome
    genes_recip: dictionary
        a dictionary of genes (CDS) and their info, with the key being the
        protein IDs and the value being a 5-element list including the
        translated sequence, the raw nucleotide sequence, the start and end
        positions in the recipient genome
    donor_genome_fp: string
        file path to donor genome
    recip_genome_fp: string
        file path to recipient genome
    output_dir: string
        path to output directory
    percentage_hgts: float
        percentage of HGT genes to simulate (of total genes in recipient
        genome)
    orthologous_rep_prob: float
        rate of orthologous replacement HGTs
    log_f: file handler
    threads: integer
        number of threads to use

    Returns
    -------        

    """
    # output dir for OrthoFinder results
    proteomes_dir = join(output_dir, "proteomes")
    if not exists(proteomes_dir):
        makedirs(proteomes_dir)

    # write donor and recipient genes to file
    if verbose:
        sys.stdout.write("Write donor and recipient genes to file.\n")
    genes_donor_fp = join(
        proteomes_dir, "%s_donor.faa" % basename(donor_genome_fp))
    genes_recip_fp = join(
        proteomes_dir, "%s_recip.faa" % basename(recip_genome_fp))
    with open(genes_donor_fp, 'w') as genes_donor_f:
        for gene in genes_donor:
            genes_donor_f.write(">%s\n%s\n" % (gene, genes_donor[gene][0]))
    with open(genes_recip_fp, 'w') as genes_recip_f:
        for gene in genes_recip:
            genes_recip_f.write(">%s\n%s\n" % (gene, genes_recip[gene][0]))

    # simulate orthologous replacement
    if orthologous_rep_prob > 0.0:
        if verbose:
            sys.stdout.write("\tSimulate orthologous replacement HGTs ...")
        launch_orthofinder(proteomes_dir, threads)
        date = time.strftime("%c").split()
        results_dir = join(
            proteomes_dir, "Results_%s%s" % (date[1], date[2]),
            "WorkingDirectory")
        species_ids, sequence_ids, orthologous_groups =\
            parse_orthofinder(results_dir)
        # no orthologs found, exit orthologous replacement simulation
        if orthologous_groups:
            seq_recip = simulate_orthologous_rep(genes_donor,
                                                 seq_donor,
                                                 genes_recip,
                                                 seq_recip,
                                                 sequence_ids,
                                                 orthologous_groups,
                                                 orthologous_rep_prob,
                                                 percentage_hgts,
                                                 log_f)
        else:
            print("WARNING: No orthologous genes found between donor and "
                  "recipient genome, continuing to novel gene acquisition "
                  "HGT simulation (if option selected).")
        if verbose:
            sys.stdout.write(" done.\n")
    # simulate novel gene acquisition
    if float(1-orthologous_rep_prob) > 0.0:
        if verbose:
            sys.stdout.write("\tSimulate novel gene acquisition HGTs ...")
        seq_recip = simulate_novel_acq(genes_donor,
                                       seq_donor,
                                       genes_recip,
                                       seq_recip,
                                       orthologous_rep_prob,
                                       percentage_hgts,
                                       log_f)
        if verbose:
            sys.stdout.write(" done.\n")
    return write_results(
        genes_donor=genes_donor,
        donor_genome_fp=donor_genome_fp,
        genes_recip=genes_recip,
        recip_genome_fp=recip_genome_fp,
        seq_donor=seq_donor,
        seq_recip=seq_recip,
        output_dir=output_dir)


def simulate_genbank(donor_genbank_fp,
                     recipient_genbank_fp,
                     output_dir,
                     percentage_hgts,
                     orthologous_rep_prob,
                     log_f,
                     threads,
                     verbose=False):
    """ Simulate HGTs using genuine genomes (GenBank format).

    Parameters
    ----------
    donor_genbank_fp: string
        file path to genome (donor of HGTs) in GenBank format
    recipient_genbank_fp: string
        file path to genome (recipient of HGTs) in GenBank format
    output_dir: string
        output directory path
    percentage_hgts: float
        percentage of HGT genes to simulate (of total genes in recipient
        genome)
    orthologous_rep_prob: float
        rate of orthologous replacement HGTs
    log_f: file handler
    threads: integer
        number of threads to use
    """
    if verbose:
        sys.stdout.write("Parsing donor GenBank record ...\n")
    seq_donor, genes_donor = extract_genbank(donor_genbank_fp, verbose)
    if verbose:
        sys.stdout.write("\tDone.\n")
    if verbose:
        sys.stdout.write("Parsing recipient GenBank record ...\t")
    seq_recip, genes_recip = extract_genbank(recipient_genbank_fp, verbose)
    if verbose:
        sys.stdout.write("\tDone.\n")

    return simulate_hgts(
        seq_donor=seq_donor,
        genes_donor=genes_donor,
        seq_recip=seq_recip,
        genes_recip=genes_recip,
        donor_genome_fp=donor_genbank_fp,
        recip_genome_fp=recipient_genbank_fp,
        output_dir=output_dir,
        percentage_hgts=percentage_hgts,
        orthologous_rep_prob=orthologous_rep_prob,
        log_f=log_f,
        threads=threads,
        verbose=verbose)


def simulate_azad_lawrence(donor_artificial_fp,
                           donor_artificial_annotation_fp,
                           recip_artificial_fp,
                           recip_artificial_annotation_fp,
                           output_dir,
                           percentage_hgts,
                           orthologous_rep_prob,
                           log_f,
                           threads,
                           verbose=False):
    """ Simulate HGTs using artificial genomes by Azad and Lawrence, 2005.

    Parameters
    ----------
    donor_artificial_fp: string
        file path to artificial genome (donor of HGTs)
    recip_artificial_fp: string
        file path to artificial genome (recipient of HGTs)
    output_dir: string
        output directory path
    percentage_hgts: float
        percentage of HGT genes to simulate (of total genes in recipient
        genome)
    orthologous_rep_prob: float
        rate of orthologous replacement HGTs
    log_f: file handler
    threads: integer
        number of threads to use
    """
    seq_donor, genes_donor = extract_azad_lawrence(
        donor_artificial_fp, donor_artificial_annotation_fp, 'donor')
    seq_recip, genes_recip = extract_azad_lawrence(
        recip_artificial_fp, recip_artificial_annotation_fp, 'recip')

    return simulate_hgts(
        seq_donor=seq_donor,
        genes_donor=genes_donor,
        seq_recip=seq_recip,
        genes_recip=genes_recip,
        donor_genome_fp=donor_artificial_fp,
        recip_genome_fp=recip_artificial_fp,
        output_dir=output_dir,
        percentage_hgts=percentage_hgts,
        orthologous_rep_prob=orthologous_rep_prob,
        log_f=log_f,
        threads=threads)


@click.command()
@click.option('--donor-genbank-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to genome (donor of HGTs) in GenBank format')
@click.option('--recipient-genbank-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to genome (recipient of HGTs) in GenBank '
                   'format')
@click.option('--donor-artificial-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to artificial genome (donor of HGTs) by Azad '
                   'and Lawrence, 2005')
@click.option('--donor-artificial-annotation-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to artificial genome annotation (donor of '
                   'HGTs) by Azad and Lawrence, 2005')
@click.option('--recipient-artificial-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to artificial genome (recipient of HGTs) by '
                   'Azad and Lawrence, 2005')
@click.option('--recipient-artificial-annotation-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='File path to artificial genome annotation (recipient of '
                   'HGTs) by Azad and Lawrence, 2005')
@click.option('--output-dir', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory path')
@click.option('--percentage-hgts', required=False, type=float, default=0.05,
              show_default=True, help='Percentage of HGT genes to simulate '
                                      '(of total genes in recipient genome)')
@click.option('--orthologous-rep-prob', required=False, type=float,
              default=0.5, show_default=True,
              help='Probability of orthologous replacement HGT (the remainder'
                ' will be HGTs in the form of novel gene acquisition)')
@click.option('--simulation-type', required=True,
              type=click.Choice(['genbank', 'azad_lawrence']),
              help='HGTs can be simulated using GenBank files (genuine '
                   'genomes) or using the annotation files by Azad & '
                   'Lawrence, 2005 (artificial genomes)')
@click.option('--threads', required=False, type=int, default=1,
              show_default=True, help='Number of threads to use')
@click.option('--verbose', required=False, type=bool, default=False,
              show_default=True, help='Run program in verbose mode')

def _main(donor_genbank_fp,
          recipient_genbank_fp,
          donor_artificial_fp,
          donor_artificial_annotation_fp,
          recipient_artificial_fp,
          recipient_artificial_annotation_fp,
          output_dir,
          percentage_hgts,
          orthologous_rep_prob,
          simulation_type,
          threads,
          verbose):
    """ Simulate HGTs by combining genes from two genomes.
    """
    if verbose:
        sys.stdout.write("Begin simulation.\n")
    # check correct input files given for corresponding simulation type 
    if simulation_type == 'genbank':
        if (donor_genbank_fp == None or
                recipient_genbank_fp == None):
            raise ValueError("The donor and recipient GenBank file paths are "
                             "required with genome-type GenBank")
    else:
        if (donor_artificial_fp is None or
                donor_artificial_annotation_fp is None or
                recipient_artificial_fp is None or
                recipient_artificial_annotation_fp is None):
            raise ValueError("The donor and recipient artificial genome file "
                             "paths and their annotation file paths are "
                             "required with genome-type azad_lawrence")

    if not exists(output_dir):
        makedirs(output_dir)
    log_fp = join(output_dir, "log.txt")
    with open(log_fp, 'w') as log_f:
        if simulation_type == 'genbank':
            if verbose:
                sys.stdout.write("Simulate HGTs in genuine genomes.\n")
            simulate_genbank(
                donor_genbank_fp=donor_genbank_fp,
                recipient_genbank_fp=recipient_genbank_fp,
                output_dir=output_dir,
                percentage_hgts=percentage_hgts,
                orthologous_rep_prob=orthologous_rep_prob,
                log_f=log_f,
                threads=threads,
                verbose=verbose)
        else:
            if verbose:
                sys.stdout.write("Simulate HGTs in artificial genomes.\n")
            simulate_azad_lawrence(
                donor_artificial_fp=donor_artificial_fp,
                donor_artificial_annotation_fp=\
                    donor_artificial_annotation_fp,
                recip_artificial_fp=recipient_artificial_fp,
                recip_artificial_annotation_fp=\
                    recipient_artificial_annotation_fp,
                output_dir=output_dir,
                percentage_hgts=percentage_hgts,
                orthologous_rep_prob=orthologous_rep_prob,
                log_f=log_f,
                threads=threads,
                verbose=verbose)


if __name__ == "__main__":
    _main()