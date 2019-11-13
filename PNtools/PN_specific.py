import pymaid
import pandas as pd
import numpy as np
import itertools
from . import misc

def PN_axon_prune(neurons,vols = None, resize = 1):
    """ Rough pruning of PNs to the axon.

    This is done by first isolating the longest neurite, and finding either the first branch point in the final neuropil which is
    inervated by the longest neurite, or the last branch point along the primary neurite if it does not end in a volume.
    The second step is to remove nodes within the right and left AL.

    Parameters
    ----------

    neurons  :   CatmaidNeuron | CatmaidNeuronList
                 Projection neuron(s) which you wish to prune to the axon

    vols     :  dict
                Dictionary of all FAFB neuropils (see PNtools.FAFB_vols function). If None (default), will run PNtools.FAFB_vols.

    resize   :  int
                Scaling factor by which to resize the Antenal lobe volumes. 1 (no scaling) by default. Useful if you want to remove
                a broader area of dendrites, but runs the risk of returning neurons with no nodes.

    Returns
    -------

    pruned   :  CatmaidNeuron | CatmaidNeuronList
                Neuron(s) pruned aproximately to leave just the axon.

    """

    if isinstance(neurons,pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    if vols is None:
        vols = misc.FAFB_vols()

    AL_R = pymaid.get_volume('AL_R_manual')
    AL_L = pymaid.get_volume('AL_L')
    pruned = pymaid.CatmaidNeuronList([])

    for N in neurons:

        pymaid.reroot_neuron(N,N.soma,inplace = True)
        ##Bit 1

        # get neurite of whole neuron
        neurite = pymaid.longest_neurite(N)
        # get end node of neurite in a volume...
        location = []
        limit = 0
        while (len(location) == 0 and limit <= 15):
            limit += 1
            last_node = list(set(neurite.nodes.treenode_id.values) - set(neurite.nodes.parent_id.values))[0]
            last_node_coords = neurite.nodes.loc[neurite.nodes.treenode_id == last_node][['x','y','z']]
            # determine which volume it is in
            location = misc.point_in_vol(last_node_coords,vols)

        if limit > 15:
            print(N.skeleton_id)
            continue
        # get branch points on neurite
        # get all branch nodes within the neuron
        dist = set(N.nodes.loc[N.nodes.type == 'branch'].treenode_id.values)
        # keep only those which are along the primary neurite
        dist = dist.intersection(set(neurite.nodes.treenode_id))
        # subset to the volume

        # coords of nodes
        # get the xyz coords of these and check which are in the final volume
        coords = pd.DataFrame()
        for i in dist:
            coords = pd.concat([coords,N.nodes.loc[N.nodes.treenode_id == i][['x','y','z']]], sort = False)

        # binary list showing which branches are in the final volume.
        keep = list(pymaid.in_volume(coords, pymaid.get_volume(location)))

        if sum(keep) > 0:
            dist = list(itertools.compress(dist,keep))


            # Find the branch closest to the root
            # of the remaining nodes, get distance to root
            dist = pd.DataFrame.from_dict({str(n) : pymaid.dist_between(N,n,N.root) for n in dist},
                                                   orient = 'index',
                                                   columns = ['Dist'])
            # Parent
            cut = neurite.nodes.loc[neurite.nodes.treenode_id.values == int(dist['Dist'].idxmin())]['parent_id']
            neurite.prune_distal_to(cut)
            # Get set of nodes along neurite between parents and root (set A)
            neurite = set(neurite.nodes.treenode_id)
        else:
                    # Find the branch closest to the root
            # of the remaining nodes, get distance to root
            dist = pd.DataFrame.from_dict({str(n) : pymaid.dist_between(N,n,N.root) for n in dist},
                                                   orient = 'index',
                                                   columns = ['Dist'])
            # Parent
            cut = neurite.nodes.loc[neurite.nodes.treenode_id.values == int(dist['Dist'].idxmax())]['parent_id']
            neurite.prune_distal_to(cut)
            # Get set of nodes along neurite between parents and root (set A)
            neurite = set(neurite.nodes.treenode_id)

        ### Bit TWO

        # Expand the AL by some factor
        AL_R.resize(resize)
        AL_L.resize(resize)
        # prune the larger AL out of the neuron
        prune = N.prune_by_volume(AL_R,mode = 'OUT', inplace = False)
        prune = prune.prune_by_volume(AL_L,mode = 'OUT', inplace = False)
        # work out which nodes to keep (set of all nodes - neurite set)
        keep = list(set(prune.nodes.treenode_id) - neurite)
        # subset neuron
        prune = pymaid.subset_neuron(prune, keep)
        pruned += prune

    return (pruned)
