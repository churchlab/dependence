"""
Recode genes.

Created: 12/4/13

Uses getk with gleb/multiplegenes branch at commit
d7eae139d8d44ef7c0712a4817172f2bb10d55b6.
https://github.com/churchlab/getk
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

from getk.aggregate_scorer import BoundedAggregateScorer
from getk.biopython_util import get_genome_record
from getk.codon_replacer import SimpleCodonReplacer
from getk.codon_usage_memex import get_ecoli_codon_usage_memex
from getk.iterators import CodonShuffler
from getk.scorer import GCScorer
from getk.scorer import SSScorer
from getk.scorer import MaintainTranslationScorer
from getk.scorer import CodonShuffleScorer
from getk.refactor_context import RefactorContext


MG1655_GENBANK= '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/mg1655/mg1655.genbank'

TYRS_ORIGINAL = 'tyrS.original.fa'
TYRS_RECODED = 'tyrS.recoded.fa'

ADK_ORIGINAL = 'adk.original.fa'
ADK_RECODED = 'adk.recoded.fa'


def extract_original_genes_from_mg1655():
    mg1655_seq_record = get_genome_record(MG1655_GENBANK)

    gene_dest_tuple_list = [
        ('tyrS', TYRS_ORIGINAL),
        ('adk', ADK_ORIGINAL),
    ]

    for gene, dest in gene_dest_tuple_list:
        for feature in mg1655_seq_record.features:
            is_gene_found = False
            if (feature.type == 'gene' and
                    feature.qualifiers['gene'][0] == gene):
                # Write the sequence.
                gene_seq_record = feature.extract(mg1655_seq_record)
                gene_seq_record.id = gene
                gene_seq_record.description = 'original'
                with open(dest, 'w') as output_fh:
                    SeqIO.write(gene_seq_record, output_fh, 'fasta')

                # Continue to next gene.
                is_gene_found = True
                break
        assert is_gene_found, '%s not found' % gene


def recode_gene(gene, original_location, recoded_dest):
    """Recode the gene and write to recoded_dest.
    """
    with open(original_location) as input_fh:
        original_seq_record = SeqIO.read(input_fh, 'fasta')

    # Add feature to the record to match SimpleCodonReplacer interface.
    feature_loc = FeatureLocation(ExactPosition(0),
            ExactPosition(len(original_seq_record)))
    original_seq_record.features = [
            SeqFeature(feature_loc, type='gene', strand=1)]

    # Build the global context object that is passed around by recoding
    # objects.
    ecoli_codon_usage_memex = get_ecoli_codon_usage_memex(randomized=True)
    params = {
        'original_seq': original_seq_record,
        'codon_usage_memex': ecoli_codon_usage_memex,
        'forbidden_codons': []
    }
    refactor_context = RefactorContext(params)

    # Iterators.
    iterators = [CodonShuffler(refactor_context)]

    # Hard scorers. Scores must be with
    hard_scorer_list = [
        GCScorer(refactor_context),
        SSScorer(refactor_context)
    ]
    hard_scorer_thresholds = [
        20, # / 100 for GC
        7, # kJ/mol
    ]

    # Put all the scorers together.
    aggregate_scorer = BoundedAggregateScorer(
            hard_scorer_list, hard_scorer_thresholds)

    # Perform the shuffling.
    codon_replacer = SimpleCodonReplacer(refactor_context, iterators,
            aggregate_scorer)
    result = codon_replacer.replace_codons_in_feature(
            original_seq_record,
            start_index=3,
            end_index=len(original_seq_record) -3)
    recoded_seq_record= result['new_seq']
    recoded_seq_record.id = gene
    recoded_seq_record.description = 'shuffled'

    # Write the result.
    with open(recoded_dest, 'w') as output_fh:
        SeqIO.write(recoded_seq_record, output_fh, 'fasta')


def analyze_recoding():
    # with open(TYRS_ORIGINAL) as fh:
    #     original = SeqIO.read(fh, 'fasta')

    # with open(TYRS_RECODED) as fh:
    #     recoded = SeqIO.read(fh, 'fasta')

    with open(ADK_ORIGINAL) as fh:
        original = SeqIO.read(fh, 'fasta')

    with open(ADK_RECODED) as fh:
        recoded = SeqIO.read(fh, 'fasta')

    # for i in range(0, len(original), 3):
    #     if str(original[i:i+3].seq) == str(recoded[i:i+3].seq):
    #         print i
    #         print str(original[i:i+3].seq)

    match = 0.0
    for i in range(0, len(original)):
        if str(original.seq[i]) == str(recoded.seq[i]):
            match += 1
    print match / len(original)


def main():
    gene_orig_recoded_tuple_list = [
        ('tyrS', TYRS_ORIGINAL, TYRS_RECODED),
        ('adk', ADK_ORIGINAL, ADK_RECODED)
    ]

    for gene, orig, recoded in gene_orig_recoded_tuple_list:
        recode_gene(gene, orig, recoded)

    analyze_recoding()


if __name__ == '__main__':
    main()
