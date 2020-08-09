from skhep.math import LorentzVector
from particle.pdgid import has_bottom, is_hadron, charge

from copy import copy

import numpy as np

DTYPE_NP = np.dtype([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('pT', 'f8'),
                     ('mass', 'f8'), ('rap', 'f8'), ('eta', 'f8'), ('theta', 'f8'),
                     ('phi', 'f8'), ('prodx', 'f8'), ('prody', 'f8'), ('prodz', 'f8'),
                     ('prodt', 'f8'), ('pdgid', 'i4'), ('status', 'i4')])

def isBottomHadron(pdgid):
    return (has_bottom(pdgid) and is_hadron(pdgid))

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
        if charge(part['pdgid']) == 0: continue
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
