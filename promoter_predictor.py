from pathlib import Path
import re
import subprocess
import shlex
from Bio import SeqIO
import random
import pandas as pd
from subprocess import CalledProcessError


def run_promoter_prediction(fastaPaths, outputPath, verbose=1):
    """
    Runs Vero's predictor on a list of fasta files.
    Note: the return positions are in 0-based index.
    
    Several limitations:
    - output folder is relative to the script's folder. Thus, we copy and clean the output files
      at the end to the correct absolute output folder.
    - Each fasta file can only contain 1 sequence. TODO: if one fasta file contains multiple
    sequences
    """

    promoterScriptFolder = '/users/lserrano/mweber/Research_cloud/promoter_predictor_Vero'
    # a limitation of Vero's script is that you can only output in the relative local folder in the script's folder.
    outputLocalFolder = '{:d}'.format(random.randint(0, 1e12))
    dfs = []

    for fastaPath in fastaPaths:
        bio = next(SeqIO.parse(fastaPath, format='fasta'))
        id1 = bio.id
        cmd = './predictPromoters.sh {} {} Promoter_PribnowMod.txt matrix_motif_35.txt EnergyTable2.txt'\
              .format(outputLocalFolder, str(fastaPath))
        if verbose >= 1:
            print(cmd)
        cmd = shlex.split(cmd)
        try:
            cmd_output = subprocess.check_output(cmd, cwd=promoterScriptFolder)
            cmd_output = re.sub(r'\\n','\n',str(cmd_output))
        except CalledProcessError:
            print("ERROR")
        if verbose >= 1:
            print(cmd_output)

        outputResLocalFolder = '{}_results'.format(outputLocalFolder)
        promoterOutputPath = Path(promoterScriptFolder) / outputResLocalFolder
        outputFiles = list(promoterOutputPath.glob('*'))
        if verbose >= 2:
            print("outputFiles:\n", outputFiles)

        df = (pd.read_table(promoterOutputPath / '{}_final_plus.txt'.format(outputLocalFolder))
              .rename(columns={'result_p[, 2]':'score'}))
        df2 = (pd.read_table(promoterOutputPath / '{}_final_minus.txt'.format(outputLocalFolder))
               .rename(columns={'result_m[, 2]':'score'}))
        df = pd.concat([df, df2])
        df['id'] = id1
        df = df.rename(columns={
            'Pos':'pos',
            'Score_-45':'score_-45',
            'Score_Energy':'score_energy',
            'Base_penalties':'base_penalties',
            'Score_Pribnow':'score_Pribnow',
            'Strand':'strand',
            'Score_-35':'score_-35',
            'Dist_-35':'dist_-35'})
        # Convert to 0-based index
        df['pos'] = df['pos'] - 1
        df.to_csv(outputPath / '{}_promoter_prediction.csv'.format(id1))
        dfs.append(df)

        # clean up files in the script folder
        [p.unlink() for p in promoterOutputPath.glob('*')]
        promoterOutputPath.rmdir()
        
    return pd.concat(dfs)
