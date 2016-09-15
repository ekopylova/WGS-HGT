#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# FASTA to Phylip converter

import click

from skbio import TabularMSA, Protein


@click.command()
@click.option('--input-fasta-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help="Input FASTA alignment file")
@click.option('--output-phylip-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help="Output Phylip alignment file")
def main(input_fasta_fp,
         output_phylip_fp):
    """Convert FASTA to Phylip alignment
    """
    msa = TabularMSA.read(input_fasta_fp, constructor=Protein)
    msa.reassign_index(minter='id')
    msa.write(output_phylip_fp, format='phylip')


if __name__ == "__main__":
    main()