#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Bin genomes by taxonomies

import sys
from os.path import exists, join, basename
from os import makedirs
import click
from collections import OrderedDict
import bz2


def get_genome_paths(scores_average_fp, scores_repophlan_fp):
    """Return taxonomy strings and their counts.

    Parameters
    ----------
    scores_average_fp: string
        Filepath to average quality scores computed for all quality scores
        reported by (Land M.L. et al., 2014) (output of
        https://github.com/tanaes/script_bin/blob/master/filter_repophlan.py)
    scores_repophlan_fp: string
        Filepath to quality scored genomes, algorithm described in
        (Land M.L. et al., 2014) and implementation in
        https://bitbucket.org/nsegata/repophlan/src/7ebcbea05bc8baacaa988e\
        7ab5e14c8206250f60/screen.py?at=default&fileviewer=file-view-default)

    Notes
    -----
    The file scores_average_fp must contain only genomes that passed the
    quality filter.
    """
    genomes_paths = {}
    # Get genome IDs (passing quality filter)
    with open(scores_average_fp) as scores_average_f:
        next(scores_average_f)
        for line in scores_average_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            faa_lname = line[1]
            if genome_id not in genomes_paths:
                genomes_paths[genome_id] = faa_lname
            else:
                raise ValueError("Duplicate genome IDs %s" % genome_id)
    # Get taxonomies based on used genome IDs
    taxonomies = {}
    genomes_without_taxonomy = []
    with open(scores_repophlan_fp) as scores_repophlan_f:
        # header
        line = scores_repophlan_f.readline()
        line = line.strip().split('\t')
        tax_idx = line.index('taxonomy')
        for line in scores_repophlan_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            # only want tax_ids for genomes passing quality filter
            if genome_id in genomes_paths:
                taxonomy = line[tax_idx]
                if (('k__Bacteria' not in taxonomy) and
                        ('k__Archaea' not in taxonomy) and
                        ('k__Viruses' not in taxonomy) and
                        ('k__Viroids' not in taxonomy)):
                    genomes_without_taxonomy.append(genome_id)
                if taxonomy not in taxonomies:
                    taxonomies[taxonomy] = [genome_id]
                else:
                    taxonomies[taxonomy].append(genome_id)
    return taxonomies, genomes_paths, genomes_without_taxonomy


def copy_files(taxonomy,
               genome_ids,
               genome_paths,
               max_genomes,
               working_dp):
    """Copy all genomes in a single taxonomy into a new directory.

    Parameters
    ----------
    taxonomy: string
        Taxonomy string under which all genomes fall
    genome_ids: set
        All genome ids with same taxonomy
    genome_paths: dictionary    
        All genome ids (keys) and their paths (values) passing 32,000
        quality genomes filter
    max_genomes: integer
        Maximum number of genomes per species to include in tree
    working_dp: string
        Working directory

    """
    species = taxonomy.split('|')[-1]
    if 's__' not in species:
        raise ValueError("Missing species level: %s" % taxonomy)
    for i, genome_id in enumerate(genome_ids):
        # a path exists for genome in bin
        if genome_id in genome_paths:
            genome_fp = genome_paths[genome_id]
            zipfile = bz2.BZ2File(genome_fp)
            data = zipfile.read()
            genome_fp_no_bz2 = join(working_dp, basename(genome_fp[:-4]))
            print("[subtrees_compute] %s. writing %s " % (i, genome_fp_no_bz2))
            with open(genome_fp_no_bz2, 'wb') as genome_fp_no_bz2_f:
                genome_fp_no_bz2_f.write(data)
        else:
            print("[subtrees_compute] %s. skipping %s " % (i, genome_id))    


@click.command()
@click.option('--repophlan-scores-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help="Output of RepoPhlAn's repophlan_get_microbes.py")
@click.option('--repophlan-scores-average-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help="Filepath to average quality scores computed for "
                   "all quality scores reported by (Land M.L. et al., 2014)")
@click.option('--working-dp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help="Working directory path")
@click.option('--max-genomes', required=True, type=int,
              help="Maximum number of genomes per species")
def main(repophlan_scores_fp,
         repophlan_scores_average_fp,
         working_dp,
         max_genomes):
    """Compute distant matrix for a subset of genomes from same species.
    """
    if not exists(working_dp):
        makedirs(working_dp)

    # 1) Bin genomes by taxonomies
    taxonomies, genomes_paths, genomes_without_taxonomy = get_genome_paths(
        repophlan_scores_average_fp, repophlan_scores_fp)
    if len(genomes_without_taxonomy) != 0:
        print("Genomes without taxonomy: %s" % len(genomes_without_taxonomy))

    # 2) Run MSA (MAFFT) on all bins having >max_genomes genomes
    # sort taxonomies by ascending number of genomes
    taxonomies = OrderedDict(sorted(taxonomies.items(), key=lambda x: len(x[1])))
    for key, genome_ids in taxonomies.items():
        if len(genome_ids) > max_genomes:
            print("[subtrees_compute] %s has %s genomes" % (key, len(genome_ids)))
            # Create directory with all genomes under same taxonomy
            output_dp = join(working_dp, key.replace('|','.'))
            if not exists(output_dp):
                makedirs(output_dp)
            copy_files(key, genome_ids, genomes_paths, max_genomes, output_dp)
    

if __name__ == "__main__":
    main()
