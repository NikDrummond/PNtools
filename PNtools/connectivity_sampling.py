# Upstream sheet generation function
import pymaid
import pandas as pd
import fafbseg
from . import misc

def upstream_node_check(neurons,volume = None):
    """ Check if upstream connectors have a node attached to them. If not, creates a DataFrame with URLs to v14 at site of the connector.

    Parameters
    ----------

    neurons:    CatmaidNeuron | CatmaidNeuronList
                A neuron or neuron list to check

    volume:     Volume
                A volume to prune the neuron(s) to.

    Returns
    -------
    DataFrame
                Data frame for each connector ID which does not have an upstream node, and URL to the connectors possiton in v14.
                If all connectors have an upstream node, will return None.
    """

    # if volume provided, prune neuron
    if volume is not None:
        neurons = pymaid.in_volume(neurons,volume)

    if isinstance(neurons,pymaid.CatmaidNeuron):
        neurons = pymaid.CatmaidNeuronList(neurons)

    # get missing bits
    missing_pre = pd.DataFrame()
    for i in neurons:
        conn = pymaid.get_connector_details(i.postsynapses)
        missing = conn[conn.presynaptic_to_node.isnull()].connector_id.values
        missing = i.connectors[i.connectors.connector_id.isin(missing)]
        missing_pre = pd.concat([missing_pre,missing])
    missing_pre = missing_pre[['connector_id','x','y','z']]

    # generate URLs
    coords = [[missing_pre.loc[i].x , missing_pre.loc[i].y, missing_pre.loc[i].z] for i in missing_pre.index]
    missing_pre['Manual_URL'] = [pymaid.url_to_coordinates(coords[i],5) for i in range(len(coords))]

    if len(missing_pre) == 0:
        return(None)
    else:
        return(missing_pre)

def upstream_sheet(neuron,volume = None,order='manual',auto_version = 'v3'):
    """ Generate a sheet with urls for all upstream neurons of a neuron.

    By default will try to order the output by number of inputs. If you have connectors without an upstream
    node, the function will give you the option to retun a DataFrame with URLs to those connectors
    instead.

    IMPORTANT: if using order = 'auto' you NEED to use fafbseg.use_google_storage, fafbseg.use_brainmaps or fafbseg.use_remote_service
    to set the way you want to fetch segmentation IDs.

    Perameters
    ----------

    neuron:         CatmaidNeuron
                    A neuron to generate the upstream sheet for

    volume:         Volume
                    A volume to prune the neuron to before generating the sheet

    order:          str
                    If 'auto' (default) will try to use the Google autoseg to order the output by the segments
                    with the greatest number of upstream points in.
                    If 'manual' will order output based on the number of inputs from the set of upstream
                    skeletons in v14.
                    if 'random' will return a randomised DataFrame.

    auto_version:   str
                    The autoseg version to use, default is 'v3' but 'v2' and 'v1' are also accepted

    Returns
    -------

    DataFrame
                    Upstream sheet of your neuron, with URLs to upstream nodes, or a sheet of inputs to
                    your neuron that do not have an upstream node.
    """

    # if volume provided, prune neuron
    if volume is not None:
        neuron = pymaid.in_volume(neuron,volume)

    # Get connectors
    conn = pymaid.get_connector_details(neuron.postsynapses)
    # get upstream nodes/neurons
    upstream = pymaid.find_treenodes(treenode_ids = list(conn.presynaptic_to_node.values))
    upstream['connector_id'] = [conn.loc[((conn.presynaptic_to_node == upstream.loc[i].treenode_id)
                                      & (conn.presynaptic_to == upstream.loc[i].skeleton_id))].connector_id.values[0] for i in upstream.index]
    upstream = upstream[['skeleton_id','treenode_id','connector_id','parent_id','x','y','z']]

    if auto_version == 'v3':
        auto_version = 'v14-seg-li-190805.0'
    elif auto_version == 'v2':
        auto_version = 'v14seg-Li-190411.0'
    elif auto_version == 'v1':
        auto_version = 'v14-seg'

    # Check for connectors with no upstream node
    missing_pre = upstream_node_check(neuron,volume)

    ans = 'z'

    while ans not in ['y','n']:
        if missing_pre is not None:
            print('You have ' +  str(len(missing_pre)) + ' upstream connectors with no upstream node!')
            print('Would you rather generate a sheet with URLs to the connectors with no ustream node?')
            ans = input(prompt = '[y|n]')
        else:
            print('All incoming connectors have an upstream node, pat yourself on the back! Or not. What do I care? No one listens to me anyway')
            ans = 'n'

    # Generate bare bones output
    if ans == 'y':
        data = missing_pre[['connector_id','x','y','z']]
    else:
        data = upstream

    if ans == 'n':
        # Add list of x,y,z coords and URLs
        coords = [[data.loc[i].x , data.loc[i].y, data.loc[i].z] for i in data.index]
        data['Manual_URL'] = [pymaid.url_to_coordinates(coords[i],5) for i in range(len(coords))]
        data['Auto_URL'] = [data.Manual_URL[i].replace('v14', auto_version) for i in range(len(data.Manual_URL))]
        if order == 'auto':
            # add fragment id column
            data['Fragment_id'] = fafbseg.segmentation.get_seg_ids(coords)
            # order
            grouper = data.groupby('Fragment_id')
            N_dict = grouper.treenode_id.count().to_dict()
            N_dict[0] = 0
            data['Hits'] = data.Fragment_id.map(N_dict)
            data.sort_values('Hits', ascending=False, inplace=True)
            data.reset_index(drop = True,inplace = True)
        elif order == 'manual':
            grouper = data.groupby('skeleton_id')
            N_dict = grouper.treenode_id.count().to_dict()
            data['Hits'] = data.skeleton_id.map(N_dict)
            data.sort_values('Hits', ascending=False, inplace=True)
            data.reset_index(drop = True,inplace = True)
        elif order == 'random':
            data = data.sample(frac=1).reset_index(drop = True)
    return(data)

def connectors_in_vol(source, volumes = None, direction = 'Both', count = False):
    """ Returns the volume(s) a neuron(s) synapses are located in.

    Parameters
    ----------

    source:     CatmaidNeuron | CatmaidNeuronList | DataFrame | list
                The neuron(s) or set of connectors (with information) you wish to get volume IDs for.
                If DataFrame give, data frame MUST have columns called ['connector_id','skeleton_id','x','y','z']
                with the relevant information in.
                A list of connector IDs can be passed, in which case skids will be determined as the presynaptic neuron.

    volumes:    Volume
                The volume(s) to count connectos in. If not given, falls back to a list of all volumes used to
                make up the FAFB full neuropil volume. See `PNtools.FAFB_vols`.

    direction:  str
                'Both' (default) returns all connectors, 'Presynaptic' returns output sites, and 'Postsynaptic'
                returns input sites. This is only applicable if a neuron/neuron list is provided as input.

    count:      Bool
                If False (default) DataFrame with each connector as a row is returned, with a column for skid, and
                a column for the volume the connector is in.
                If True a DataFrame is returned with a row for each volume, and each column as a neuron, where values
                are a count of the number of connectors that neuron has within the volume.

    Returns
    -------

    DataFrame
                Either a connectors by [Skids, Volume] data frame, or a Volume by neurons data frame with connector counts
                (see count option above).

    """

    # if source is neuron/neuron list get connectors data frame
    if (isinstance(source, pymaid.CatmaidNeuron)|(isinstance(source, pymaid.CatmaidNeuronList))):
        if isinstance(source,pymaid.CatmaidNeuron):
            source = pymaid.CatmaidNeuronList(source)
        # if direction specified, subset accordingly
        if direction == 'Both':
            source = source.connectors
        elif direction == 'Presynaptic':
            source = source.presynapses
        elif direction == 'Postsynaptic':
            source = source.postsynapses

    if isinstance(source, list):

        # get node locations, order by node_id
        source = pymaid.get_node_location(source)
        source = source.sort_values(by = 'node_id', axis = 0)
        # get connector deets, order by connector id
        conn = pymaid.get_connector_details(source.node_id)
        conn = conn.sort_values(by = 'connector_id', axis = 0)
        # add skeleton id col. from connector deets to node locations etc...
        source['skeleton_id'] = conn.presynaptic_to
        source = source.rename(columns = {'node_id':'connector_id'})
        source = source[['connector_id','skeleton_id','x','y','z']]

    if volumes is None:
        volumes = misc.FAFB_vols()
    elif isinstance(volumes, pymaid.Volume):
        volumes = {volumes.name : volumes}

    skids = set(source.skeleton_id)
    skids = [str(s) for s in skids]
    # for each unique skid in source:
    data = pd.DataFrame()
    if count:
        for s in skids:
            # subset source to that skid, then get dictionary
            dictionary = pymaid.in_volume(source.loc[source.skeleton_id == s],volumes)
            data = pd.concat([data,pd.DataFrame.from_dict({n:sum(dictionary[n]) for n in dictionary.keys()},
                                          orient = 'index',
                                          columns = [str(s)])],axis = 1)
    else:
        dictionary = pymaid.in_volume(source,volumes)
        data = pd.DataFrame(data = source.skeleton_id.values,
                               index = source.connector_id,
                               columns = ['skeleton_id'])
        data['Volume'] = ""
        for k in dictionary.keys():
            indicies = [i for i, x in enumerate(dictionary[k]) if x]
            if len(indicies) > 0:
                data.iloc[indicies,[1]] = k

    return (data)
