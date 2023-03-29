import numpy as np
from Bio.Seq import Seq
from scipy.signal import argrelextrema
import RNA



def compute_hybridization_energy_along_seq(motif, seq, verbose=1, reverse_complement_motif=True):
    motifBio = Seq(motif)
    if reverse_complement_motif:
        motifSeq = str(motifBio.reverse_complement())
    else:
        motifSeq = str(motifBio)

    if verbose >= 2: print("motifSeq", motifSeq)
    l = len(motifSeq)
    if len(seq) < l:
        subSeqList = [seq]
    else:
        subSeqList = [seq[i:i+l] for i in range(len(seq) - l + 1)]
    if verbose >= 2: print("subSeqList", subSeqList)
    energyFoldList = [RNA.cofold(motifSeq + '&' + subSeq) for subSeq in subSeqList]
    return energyFoldList


def find_motif_positions_with_low_energy(seq, motif, verbose=1):
    
    energyFoldList = compute_hybridization_energy_along_seq(motif, seq, verbose=verbose)
    energyArr = np.array([e for fold, e in energyFoldList])
    energyNonZeroPos = [(i, x) for i, x in enumerate(energyArr) if x < 0]
    energyNonZero = [x for i, x in energyNonZeroPos]
    energyNonZeroIndices = [i for i, x in energyNonZeroPos]
    energyLocalMinima = argrelextrema(np.array(energyNonZero), np.less_equal)[0]
    energyMinimaIndices = [energyNonZeroIndices[i] for i in list(energyLocalMinima)]
    energyMinima = [(i, energyArr[i]) for i in energyMinimaIndices]
    return energyMinima