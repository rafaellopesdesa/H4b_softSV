from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID

from pyjet import cluster, DTYPE_EP
from pyjet.testdata import get_event

from skhep.math import LorentzVector

import numpy as np
from numpy.lib.recfunctions import append_fields, merge_arrays
import pandas as pd

from sys import argv

from utils import *

    
# We will follow ATLAS recommendations here
# No neutrinos and no muons
selection_stable = ((STATUS == 1) & ~HAS_END_VERTEX &
                    (ABS_PDG_ID != 12) & (ABS_PDG_ID != 14) & (ABS_PDG_ID != 16) & (ABS_PDG_ID != 13))


# data structure for saving
df_columns = ['nb' , 'nsv', 'hm']
category_data = []

# this can be the same
pythia = Pythia(params={'Random:seed':  '%s' % argv[2]}, config='./decay_findB.cmnd', random_state=1)
totalEvents = int(argv[1])

# but the analysis has to be a bit different
# in particular, I am interested in doing labeling like ATLAS does
for ievt, event in enumerate(pythia(events=totalEvents)):

    if ievt%100 == 0:
        print '%d/%d' % (ievt, totalEvents)
    
    array_all = event.all(return_hepmc = True)
    array_copy  = event.all()
    
    secvtx = []
    all_b = []
    for part in array_all:
        if part.pid == 35 and is_last_copy(part):
            for a_boson in find_a_boson(part):
                b_hadrons = []
                find_b_hadron(a_boson, b_hadrons)
                all_b += b_hadrons
                for b_hadron in b_hadrons:
                    sv = find_sv_hadron(b_hadron)
                    if sv is not None:
                        secvtx.append(sv)

    not_b = []
    for i,b in enumerate(array_copy):
        if deltaR(min(all_b, key = lambda k: deltaR(k, b)),b) > 1e-18:
            not_b.append(i)
    b_ghosts = np.delete(array_copy, not_b)
    ghostify(b_ghosts,1e-18)
    b_ghosts = append_fields(b_ghosts, 'bid', data=np.ones_like(b_ghosts['pdgid']), dtypes='<i4')
    
    array_st = event.all(selection_stable)
    array_st = append_fields(array_st, 'bid', data=np.zeros_like(array_st['pdgid']), dtypes='<i4')

    sequence = cluster(np.hstack([array_st,b_ghosts]), R=0.4, p=-1, ep=True)
    jets = sequence.inclusive_jets()
    selected_jets = np.array([jet for jet in sequence.inclusive_jets() if (jet.pt > 15 and np.abs(jet.eta) < 2.5)])

    bids = {}
    for jet in selected_jets:
        bids[jet] = sum([const.bid for const in jet])

    selected_secvtx = []    
    if len(selected_jets) > 0:
        for sv in secvtx:
            dr = deltaR(min(selected_jets, key = lambda k: deltaR(k,sv)),sv)
            if dr > 0.4 and sv['pT'] > 3 and sv['mass'] > 0.6:
                selected_secvtx.append(sv)

    #secondary vertices merge sometimes
    merged_secvtx = np.array(selected_secvtx,dtype=DTYPE_NP)
    merged_secvtx[merged_secvtx['pT'].argsort()]
    merged_secvtx = merge_sv(merged_secvtx, 0.75)
    
    not_b = []
    for i,jet in enumerate(selected_jets):
        if bids[jet] == 0:
            not_b.append(i)
    selected_bjets = np.delete(selected_jets, not_b)
    
    higgs_px = sum([jet.px for jet in selected_bjets]) + sum([sv['px'] for sv in merged_secvtx])
    higgs_py = sum([jet.py for jet in selected_bjets]) + sum([sv['py'] for sv in merged_secvtx])
    higgs_pz = sum([jet.pz for jet in selected_bjets]) + sum([sv['pz'] for sv in merged_secvtx])
    higgs_e = sum([jet.e for jet in selected_bjets]) + sum([sv['E'] for sv in merged_secvtx])
    higgs = LorentzVector(higgs_px,higgs_py,higgs_pz,higgs_e)
    category_data.append([len(selected_bjets), len(merged_secvtx), higgs.mass])

track_df = pd.DataFrame(category_data, columns=df_columns)
track_df.to_hdf('data_secvtx_%d.h5' % int(argv[2]), key='cats', mode='w')

print 'Saving ', len(category_data), ' jets in data_secvtx.h5'
