from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID

from pyjet import cluster, DTYPE_EP
from pyjet.testdata import get_event

# yeah, you need to do:
# pip install numpythia --user
# pip install pyjet --user

import numpy as np
from numpy.lib.recfunctions import append_fields, merge_arrays
import pandas as pd
# ok, you need numpy and pandas

from skhep.math import LorentzVector

from copy import copy
from sys import argv, stdout

from itertools import combinations, permutations

from PIDUtils import *

DTYPE_NP = np.dtype([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('pT', 'f8'),
                     ('mass', 'f8'), ('rap', 'f8'), ('eta', 'f8'), ('theta', 'f8'),
                     ('phi', 'f8'), ('prodx', 'f8'), ('prody', 'f8'), ('prodz', 'f8'),
                     ('prodt', 'f8'), ('pdgid', 'i4'), ('status', 'i4')])

def dphi(dphi):    
    delta_phi = copy(dphi)
    while (delta_phi >= np.pi): delta_phi -= 2*np.pi
    while (delta_phi < -np.pi): delta_phi += 2*np.pi
    return delta_phi
    
def deltaR(t1, t2):
    delta_eta = 0
    delta_phi = 0
    try:
        delta_eta = t1.eta - t2.eta
    except:
        try:            
            delta_eta = t1.eta - t2['eta']
        except:
            delta_eta = t1['eta'] - t2['eta']
            
    try:
        delta_phi = t1.phi - t2.phi
    except:
        try:
            delta_phi = t1.phi - t2['phi']
        except:
            delta_phi = t1['phi'] - t2['phi']
            
        
    delta_phi = dphi(delta_phi)
    return np.sqrt(delta_eta*delta_eta + delta_phi*delta_phi)

def is_in_gen_list(part, partlist):
    eps = 1e-18
    for part2 in partlist:
        if deltaR(part, part2) < eps:
            return True
    return False

def is_last_copy(part):
    for child in part.children(return_hepmc = True):
        if child.pid == part.pid:
            return False
    return True

def find_last_copy(part):
    for child in part.children(return_hepmc = True):
        if child.pid == part.pid:
            return find_last_copy(child)
    return part

def find_a_boson(part):
    for child in part.children(return_hepmc = True):
        if child.pid == 36:
            return [find_last_copy(a_boson) for a_boson in part.children(return_hepmc = True)]
        else:
            return find_a_boson(child)
    return

def find_b_hadron(part, b_hadrons):
    for child in part.children(return_hepmc = True):
        if isBottomHadron(child.pid):
            last_b = find_last_copy(child)
            if not is_in_gen_list(last_b, b_hadrons):
                b_hadrons.append(last_b)
        else:
            find_b_hadron(child, b_hadrons)
    return b_hadrons

def find_ewk_decay(b_hadron):
    for child in b_hadron.children(return_hepmc = True):
        hasBsibling = False
        for part in b_hadron.children(return_hepmc = True):
            if isBottomHadron(part.pid):
                hasBsibling = True
                break
        if not hasBsibling:
            return b_hadron.children(), b_hadron
        else:
            if isBottomHadron(child.pid):
                return find_ewk_decay(child)
    return

def find_sv_hadron(b_hadron):
    first_c_hadron, last_b_hadron = find_ewk_decay(b_hadron)
    sv_x = first_c_hadron[0]['prodx']     
    sv_y = first_c_hadron[0]['prody']
    sv_z = first_c_hadron[0]['prodz']
    sv_t = first_c_hadron[0]['prodt']
    d0 = np.sqrt(sv_x*sv_x + sv_y*sv_y)
    tracks = []
    for part in last_b_hadron.descendants():
        if part['status'] != 1: continue
        if not isCharged(part['pdgid']): continue
        if part['pT'] < 0.5: continue
        if abs(part['eta']) > 2.5: continue
        tracks.append(part)

    if len(tracks) == 0:
        return None
    else:
        sv = LorentzVector(sum([track['px'] for track in tracks]), 
                           sum([track['py'] for track in tracks]),
                           sum([track['pz'] for track in tracks]),
                           sum([track['E'] for track in tracks]))
        leadpt = max([track['pT'] for track in tracks])
        if leadpt > 2 and d0 > 0.05:
            return np.array([(sv.e, sv.px, sv.py, sv.pz, sv.pt, sv.mass, sv.rapidity, sv.eta, sv.theta(), sv.phi(), sv_x, sv_y, sv_z, sv_t, last_b_hadron.pid, 1)], dtype=DTYPE_NP)
        else:
            return None

def ghostify(partList, ghost_pt):
    partList['px'] *= ghost_pt/np.sqrt(partList['E']*partList['E']-partList['mass']*partList['mass'])
    partList['py'] *= ghost_pt/np.sqrt(partList['E']*partList['E']-partList['mass']*partList['mass'])
    partList['pz'] *= ghost_pt/np.sqrt(partList['E']*partList['E']-partList['mass']*partList['mass'])
    partList['E'] *= ghost_pt/partList['E']
    partList['pT'] = np.sqrt(partList['px']*partList['px']+partList['py']*partList['py'])
    partList['eta'] = np.arccosh(ghost_pt/partList['pT'])
    partList['phi'] = np.arctan2(partList['py'],partList['px'])

def merge(p1,p2):
    sv = LorentzVector(p1['px']+p2['px'], p1['py']+p2['py'], p1['pz']+p2['pz'], p1['E']+p2['E'])
    return np.array([(sv.e, sv.px, sv.py, sv.pz, sv.pt, sv.mass, sv.rapidity, sv.eta, sv.theta(), sv.phi(), p1['prodx'], p1['prody'], p1['prodz'], p1['prodt'], p1['pdgid'], 1)], dtype=DTYPE_NP)

    
def merge_sv(inputs, radius):
    if len(inputs) < 2:
        return inputs
    for sv in inputs:
        closest_sv = min([x for x in inputs if x != sv], key = lambda k: deltaR(k, sv))
        closest_dr = deltaR(sv, closest_sv)
        if closest_dr < radius:
            inputs = np.append(inputs, merge(sv, closest_sv))
            inputs = np.delete(inputs, np.where(inputs == sv))
            inputs = np.delete(inputs, np.where(inputs == closest_sv))
            inputs = merge_sv(inputs, radius)
            break
    return inputs
        
        

    
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
