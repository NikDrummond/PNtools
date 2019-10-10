# Various functions to determine the amount of cable a neuron has within a volume

import pymaid
import pandas as pd
import numpy as np
from tqdm import tqdm
from . import utils

@utils.has_remote_instance
def get_gloms(Side = 'Right', instance = None):
    """ Collects all of the Glomeruli volumes from CATMAID.

    Parameters
    ----------
    Side :      str
                Specifies whether to return glomeruli in the right, left, or both hemispheres of the FAFB volume. Takes string as input ('Right','Left','Both').
                Set to 'Right' by default. Can be set to 'FIB' to return the FIB glomeruli volumes from the local CATMAID instance. This requires that you are
                a - on the local Zoology network, and b - you have to pass the function the remote instance for the local CATMAID instance.

    instance:   CatmaidInstance
                Which remote instance to use to pull glomeruli form. If not give (default) will fall back to global instance.

    Retruns
    -------
    gloms :     dict
                A dictionary of glomeruli volumes to pymaid volumes for all glomeruli

    """

    if instance is None:
        instance = pymaid.utils._eval_remote_instance(instance)

    if Side == 'FIB':
        all_vols = pymaid.get_volume(remote_instance = instance)
        glom_names = [n for n in all_vols.name.values if n.startswith('FIB') and
                     not n.endswith('neuropil')]
    else:
        # Get a list of all volumes
        all_vols = pymaid.get_volume()
        # get a rough list of names we are interested in
        glom_names = [n for n in all_vols.name.values if n.startswith('v14') and
                     True not in [k in n for k in ['Lo','LC6', 'neuropil', 'LPC', 'LP_', 'right', '_ORNs']]]
        # Remove the specific names we don't want to. This is because the 'Either a pymaid neuron, or neuron list._new' meshes are more accurate, and we don't want the VP1 sub-volumes
        glom_names = [g for g in glom_names if g not in ['v14.VP1', 'v14.VP2', 'v14.VP3', 'v14.VP4', 'v14.VP5',
                                                         'v14.VP1m', 'v14.VP1l', 'v14.VP1d', 'v14.VC5', 'v14.VP1_L',
                                                         'v14.VP2_L', 'v14.VP3_L', 'v14.VP4_L', 'v14.VP5_L', 'v14.VP1m_L',
                                                         'v14.VP1l_L', 'v14.VP1d_L', 'v14.VC5_L']]
        # Sort left, right, or both sides
        if Side == 'Right':
            # Remove left hand side glomeruli
            glom_names = [g for g in glom_names if not g.endswith('_L')]
        elif Side == 'Left':
            # Remove right hand side Glomeruli
            glom_names = [g for g in glom_names if g.endswith('_L')]

    if instance is None:
        gloms = pymaid.get_volume(glom_names)
    else:
        gloms = pymaid.get_volume(glom_names,remote_instance = instance)

    # clean up glomeruli names
    if Side == 'FIB':
        gloms = {k.replace('FIB.',''): v for k, v in gloms.items()}
    else:
        gloms = {k.replace('v14.','').replace('_new',''): v for k, v in gloms.items()}

    return (gloms)



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
    ends:     DataFrame
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
    pruned :            CatmaidNeuron | CatmaidNeuronList
                        Pruned pymaid neuron / neuron list.

    """
    # If single neuron provided, change to neuron list
    if isinstance(neurons, pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    # resize volume
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


def vol_of_vol(volumes):
    """ Get the volume of volumes in CATMAID in nanometers cubed.

    Parameters
    ----------
    volumes:     Volume | dict
                 Either a pymaid volume, or a dictoinary of pymaid volumes.

    Returns
    -------
    DataFrame:   DataFrame
                 A data frame with the volume in nanomters cubed of each volume passed.

    """
    if isinstance(volumes, pymaid.Volume):
        volumes = {volumes.name: volumes}

    volumes = pd.DataFrame(data = [volumes[a].to_trimesh().volume for a in volumes],
                           index = volumes.keys(),
                           columns = ['Volume'])
    return (volumes)


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
    DataFrame:        DataFrame
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
        volumes = vol_of_vol(volumes)
        cable_mat = pd.DataFrame(data = cable_mat.values / volumes.T.values,
                                 index = cable_mat.index,
                                 columns = cable_mat.columns)

    return(cable_mat)
