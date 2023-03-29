# coding: utf-8
# Author: Marc Weber

"""
=========================================================================================
Paths
=========================================================================================
"""

from pathlib import Path

# Global variables solution
# Doesn't work because it is very difficult to use global variables in Python


# Class solution, works well but need to access an instance of class p
class p:

    def __init__(self):

        self.rootPath = Path('/users/lserrano/mweber')
        self.rootNoBackupPath = Path('/no_backup/lserrano/mweber')
        self.databasesPath = self.rootPath / 'Databases'

        self.homePath = Path('/users/lserrano/mweber')

        self.researchDropboxPath = self.homePath / 'Research_cloud'

        # Research Dropbox
        self.mycoplasmaDataPath = self.researchDropboxPath / 'Mycoplasma_pneumoniae_experimental_data'
        self.mpnAnnotationPath = self.mycoplasmaDataPath / 'Annotation'
        self.mpnAnnotWT2Path = self.mycoplasmaDataPath / 'Annotation_WT2'
        self.cterminalPath = self.researchDropboxPath / 'C-terminal'
        self.ctermELISAPath = self.cterminalPath / 'results' / '2017.01.30_Raul_cterminal_AA_construct_results'
        self.translationModelPath = self.researchDropboxPath / 'Translation_model'
        self.silacPath = self.researchDropboxPath / 'SILAC'
        self.translationErrorsPath = self.researchDropboxPath / 'Translation_errors'
        self.EFPPath = self.researchDropboxPath / 'EF-P_mutant'
        self.tmRNAPath = self.researchDropboxPath / 'tmRNA_mutant'
        self.KRASPath = self.researchDropboxPath / 'KRAS_codon_usage'

        # Databases
        self.refSeqPath = self.databasesPath / 'RefSeq/Release_2016_10'
        self.taxonomyPath = self.refSeqPath / 'Entrez_taxonomy'
        self.phylaPath = self.taxonomyPath / 'Phyla_genome_list'
        self.entrezGenePath = self.refSeqPath / 'Entrez_gene_mapping'
        self.OMAPath = self.databasesPath / "OMA"
        self.eggNOGPath = self.databasesPath / 'EggNOG4.5/eggnogdb.embl.de/download/eggnog_4.5'
        self.NCBI_COG_path = self.databasesPath / 'NCBI_COG'
        self.amiGOpath = self.databasesPath / 'AmiGO'
        self.uniprotPath = self.databasesPath / 'uniprot'
        self.paxdbPath = self.databasesPath / 'paxdb'

        # Sequencing
        self.rnaSeqDataPath = self.rootPath / 'RNA-seq_data'
        self.riboDataPath = self.rnaSeqDataPath / 'Ribosome_profiling_data'

        # Analysis data
        self.analysisDataPath = self.rootPath / 'Analysis_data'
        self.analysisDataNoBackupPath = self.rootNoBackupPath / 'Analysis_data'

        self.analysisTSSPath = self.analysisDataPath / 'TSS_sequencing_data'
        self.analysisCtermPath = self.analysisDataPath / 'C-terminal'
        self.analysisEggnogPath = self.analysisDataPath / 'C-terminal' / 'Eggnog'
        self.analysisAbundancePath = self.analysisDataPath / 'C-terminal' / 'Protein_abundance'
        self.analysisCtermDataPath = self.analysisDataPath / 'C-terminal' / 'Data'
        self.analysisCtermPlotsPath = self.analysisDataPath / 'C-terminal' / 'Plots'
        self.ctermELMseqPath = self.analysisCtermPath / 'DAM_screening'
        self.analysisTranslationModelPath = self.analysisDataPath / 'Translation_model2'
        self.translationModelSimulationOutputPath = \
            self.analysisDataPath / 'Translation_model2' / 'MPN_model' / 'Simulation_output'
