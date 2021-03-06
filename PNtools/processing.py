# Various functions to determine the amount of cable a neuron has within a volume

import pymaid
import pandas as pd
import numpy as np
from tqdm import tqdm
import itertools
from . import utils
from . import misc

def ends_matrix(neurons, volumes, as_mask = False):
    """ Return a count of end nodes the neuron(s) have within the given volume(s).

    Parameters
    ----------
    neurons :   CatmaidNeuron | CatmaidNeuronList
                Input neuron(s) of interest.

    volumes :   Volume | dict
                Either a single pymaid volume, or a dictionary of volumes.

    as_mask:    Bool
                If True, returns a boolian array for use as a mask. False by default.

    Returns
    -------
    DataFrame
              Returns an n by volumes pandas data frame, with the number of leaf nodes within the volume.
              Note: this ignores ends tags, and returns the count of every node with no child.

    """
    # If n is a single neuron, convert it to a neuron list
    if isinstance(neurons, pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    # Check if volumes is a single volume, and change to a dictionary if not.
    if isinstance(volumes, pymaid.Volume):
          volumes = {volumes.name : volumes}

    # Get matrix of end nodes within volumes

    end_mat = pd.DataFrame()
    for i in tqdm(neurons):
        # get count of end nodes
        open_ends = set(i.nodes.treenode_id.values) - set(i.nodes.parent_id.values)
        dictionary = pymaid.in_volume(i.nodes.loc[i.nodes['treenode_id'].isin(open_ends)][['x','y','z']], volumes)
        counts = pd.DataFrame(data = [np.sum(dictionary[x]) for x in dictionary.keys()],
                 index = dictionary.keys(),
                 columns = [i.skeleton_id])
        # concatonate into output matrix
        end_mat = pd.concat([end_mat, counts], axis = 1, sort = True)

    if as_mask == True:
        end_mat = end_mat.astype(bool)

    return end_mat.T

def path_to_root(node,neuron):
    """  Get path between a node and the neurons root """
    target = neuron.soma
    if isinstance(neuron,pymaid.CatmaidNeuron):
        neuron = neuron.graph
    path = []
    while node and node !=target:
        path.append(node)
        node = next(neuron.successors(node), None)

    return path

def first_branch(neurons, volume = None):
    """ Return node ID(s) for the parent node of the first branch point in the neuron(s)

    Parameters
    ----------

    neurons     CatmaidNeuron | CatmaidNeuronList
                A neuron, or neuron list. If neuron, returns single branch point.

    Returns
    -------

    list
                List of nodes IDs closest to neuron(s) root.


    """
    if isinstance(neurons, pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    nodes = []
    for i in neurons:

        # get the primary neurite
        neurite = pymaid.longest_neurite(i)
        # if volume provided, prune to that volume.
        if volume is not None:
            i.prune_by_volume(volume, inplace = True)
        # get all branch nodes in neuron
        dist = set(i.nodes.loc[i.nodes.type == 'branch'].treenode_id.values)
        # get the intersection, so branch points along the primary neurite
        dist = dist.intersection(set(neurite.nodes.treenode_id))
        # create a data frame with the distances to the root
        dist = pd.DataFrame.from_dict({str(n) : pymaid.dist_between(i,n,i.root) for n in dist},
                                       orient = 'index',
                                       columns = ['Dist'])
        nodes.append(i.nodes.loc[i.nodes.treenode_id.values == int(dist['Dist'].idxmin())]['parent_id'].values[0])
    return nodes

def pruning(neurons, volume, version = 'new', vol_scale = 1, prevent_fragments = False):
    """ Prunes a neuron to a volume in a manner which attempts to limit the neuron to cable which is likely to synapse.
    Parameters
    ----------
    neurons :           CatmaidNeuron | List
                        pymaid neuron or neuron list to be pruned
    volume :            Volume
                        CATMAID volume to prune the neuron to.
    version :           str
                        Defult to 'new', but can be 'old' if wished. See notes.
    vol_scale:          int
                        integer, determining how to reseize the volume mesh being passed. 1 (default) doesn't change the size,
                        0.5 would halve the size, 1.5 would increase the volume by 50% etc
    prevent_fragments:  Bool
                        If True, returns a single complete subgraph, if False (default) will potentially return a fragmented
                        neuron. The fragmented neuron will likely be better pruned, but depending on further analysis a
                        complete sub graph may be wanted.
    Returns
    -------
    CatmaidNeuron | CatmaidNeuronList
                        Pruned pymaid neuron / neuron list.
    """
    # If single neuron provided, change to neuron list
    if isinstance(neurons, pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    # resize volume if needed
    if vol_scale is not 1:
        volume.resize(vol_scale,inplace = True)
    # Initilise neuron lists for pruned inhibitory and excitatory neurons
    pruned = pymaid.CatmaidNeuronList([])
    # loop and prune
    if version == 'old':
        for i in tqdm(neurons):
            # prune the current neuron in this iteration to the AL
            i.reroot(i.soma, inplace=True)
            current = i.prune_by_volume(volume, prevent_fragments=True, inplace=False)
            # prune by strahler
            current = pymaid.prune_by_strahler(current, to_prune=slice(-1, None), inplace=False)
            # add to initialised neuron list
            pruned += i
    elif version == 'new':
        for i in tqdm(neurons):
            i.reroot(i.soma, inplace = True)
            # get the longest neurite
            neurite = pymaid.longest_neurite(i)
            # get the last node in the longest neurite
            last_node = list(set(neurite.nodes.treenode_id.values) - set(neurite.nodes.parent_id.values))
            # prune the neuron to the volume of interest as a complete graph
            vol_prune = i.prune_by_volume(volume,prevent_fragments = True, inplace = False)
            # Get set of all branch points in the neuron
            dist = set(i.nodes.loc[i.nodes.type == 'branch'].treenode_id.values)
            # Get intersection of this set with set of nodes in primary neurite and vol_prune
            dist = dist.intersection(neurite.nodes.treenode_id.values,vol_prune.nodes.treenode_id.values)
            # Get the distance of each branch node from the root
            dist = pd.DataFrame.from_dict({str(n) : pymaid.dist_between(i,n,i.soma) for n in dist},
                               orient = 'index',
                               columns = ['Dist'])
            # if the primary neurite ends in volume of interest
            if pymaid.in_volume(i.nodes.loc[i.nodes['treenode_id'] == last_node[0]][['x','y','z']],volume)[0]:
                # get the parent of the branch node closest to the root
                cut = neurite.nodes.loc[neurite.nodes.treenode_id.values == int(dist['Dist'].idxmin())]['parent_id']
                neurite.prune_distal_to(cut,inplace = True)
                # remove everything proximal to the root
                subset = list(set(vol_prune.nodes.treenode_id) - set(neurite.nodes.treenode_id))
                pymaid.subset_neuron(vol_prune,subset,clear_temp = True,inplace = True, prevent_fragments = prevent_fragments)
            # if the neurite does not end in the volume, remove primary neurite
            else:
                # subtract the entire primary neurite
                # subtract the neurite nodes from the cut neuron nodes:
                subset = list(set(vol_prune.nodes.treenode_id) - set(neurite.nodes.treenode_id))
                pymaid.subset_neuron(vol_prune,subset,clear_temp = True,inplace = True, prevent_fragments = prevent_fragments)
            pruned += vol_prune

    return (pruned)

def cable_length_matrix(neurons, volumes, mask = None, Normalisation=None):
    """ Matrix of neuron cable length (nanometers) within volume(s)

    If mask is provided the retruned matrix will return 0 cable where the mask is False.

    Parameters
    ----------
    neurons:          CatmaidNeuron | CatmaidNeuronList
                      A pymaid neuron or neuron list

    volumes:          Volume | dict
                      A single pymaid volume, or a dictionary of volumes.

    mask:             Boolian Array
                      A neuron(s) by volume(s) boolian array to use as a mask.
                      Defaults to set volumes with neurons with cable in but no end nodes to 0 cable.
                      NOTE - This is applied BEFORE normalisation

    Normalisation:    str
                      Normaliseation method to use. If not set, returns raw cable length.
                      If 'Neuron' normalises cable length by sum of all cable in all passed volumes
                      If 'Volume' normalises cable length by the volume of the volumes passed

    Returns
    -------
    DataFrame
                      Returns a Neuron(s) by Volume(s) data frame with the length of neurons cable within each volume.

    """

    # Coerce stuff to how we want if single neuron or volume given
    if isinstance(neurons, pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)
    if isinstance(volumes, pymaid.Volume):
        volumes = {volumes.name: volumes}

    # cut neuron list to within volumes
    res = pymaid.in_volume(neurons,volumes)
    # Create a data frame of cable lengths
    cable_mat = pd.DataFrame.from_dict({g: {i.skeleton_id: i.cable_length for i in res[g]} for g in res})

    # Masking
    if mask is not None:
        cable_mat = cable_mat.where(mask).fillna(0)

    # Normaliseation
    if Normalisation == 'Neuron':
        cable_mat = (cable_mat.T / cable_mat.sum(axis=1)).T
    elif Normalisation == 'Volume':
        volumes = misc.vol_of_vol(volumes)
        cable_mat = pd.DataFrame(data = cable_mat.values / volumes.T.values,
                                 index = cable_mat.index,
                                 columns = cable_mat.columns)

    return(cable_mat)
