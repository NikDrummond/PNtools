import pymaid
import pandas as pd
import numpy as np
from . import utils

def vol_of_vol(volumes):
    """ Get the volume of volumes in CATMAID in nanometers cubed.

    Parameters
    ----------
    volumes:     Volume | dict
                 Either a pymaid volume, or a dictoinary of pymaid volumes.

    Returns
    -------
    DataFrame
                 A data frame with the volume in nanomters cubed of each volume passed.

    """
    if isinstance(volumes, pymaid.Volume):
        volumes = {volumes.name: volumes}

    volumes = pd.DataFrame(data = [volumes[a].to_trimesh().volume for a in volumes],
                           index = volumes.keys(),
                           columns = ['Volume'])
    return (volumes)

@utils.has_remote_instance
def FAFB_vols(list = False):
    """ Returns core neuropils used to generate the FAFB mesh.

    Based on the `FAFBNP.surf$RegionList` function in the elmr R package.

    Parameters
    ----------

    list:  Bool
            (optional) If True, function prints a list of volumes, rather than fetch them. False by default.

    Returns
    -------

    dict | list
            By default, returns a dictionary of the "core" FAFB volumes used to create the full neuropil mesh.
            if list == True, just retruns a list of "core" neuropil names.

    """
    # Can't think of a smart way to do this, so pulled names of the volumes from the elmr package in R and just have them as a list
    eugh = ["AME_R","LO_R","NO","BU_R","PB","LH_R","LAL_R","SAD","CAN_R","AMMC_R","ICL_R",
            "VES_R","IB_R","ATL_R","CRE_R","MB_PED_R","MB_VL_R","MB_ML_R","FLA_R","LOP_R",
            "EB","AL_R","ME_R","FB","SLP_R","SIP_R","SMP_R","AVLP_R","PVLP_R","WED_R","PLP_R",
            "AOTU_R","GOR_R","MB_CA_R","SPS_R","IPS_R","SCL_R","EPA_R","GNG","PRW","GA_R",
            "AME_L","LO_L","BU_L","LH_L","LAL_L","CAN_L","AMMC_L","ICL_L","VES_L","IB_L",
            "ATL_L","CRE_L","MB_PED_L","MB_VL_L","MB_ML_L","FLA_L","LOP_L","AL_L","ME_L",
            "SLP_L","SIP_L","SMP_L","AVLP_L","PVLP_L","WED_L","PLP_L","AOTU_L","GOR_L",
            "MB_CA_L","SPS_L","IPS_L","SCL_L","EPA_L","GA_L"]

    if list:
        return (eugh)
    else:
        vols = pymaid.get_volume(eugh)
        return (vols)

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
    dict
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
                                                         'v14.VP1_new', 'v14.VP1_L', 'v14.VP2_L', 'v14.VP3_L',
                                                         'v14.VP4_L', 'v14.VP5_L', 'v14.VP1m_L', 'v14.VP1l_L',
                                                         'v14.VP1d_L']]
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

def calc_ltk(x):
    """ Calculate Lifetime Kurtosis from a list."""

    # Make sure there are no NaNs:
    x = x[ ~np.isnan(x)]
    # ltk:
    ltk = sum( ( ( x-x.mean() ) / x.std() ) **4) / len(x) - 3

    return ltk


def calc_lts(x):
    """ Calculate Lifetime Sparseness from a list."""

    # Make sure there are no NaNs:
    x = x[ ~np.isnan(x)]

    N = len(x)
    # lts:
    lts = 1 / (1 - (1 / N) ) * ( 1 - ( sum( x / N ) **2 / sum( x**2 / N ) ) )

    return lts
