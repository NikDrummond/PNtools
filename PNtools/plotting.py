# Wrapper functions for useful plotting
import matplotlib.pyplot as plt
import pandas as pd
from math import pi

def radar_plot(data, subsets = None, size = (10,10)):
    """ Create a radar, or spider plot to show, for example, glomeruli representations in volumes

    Parameters
    ----------

    data:       DataFrame
                Columns are 'spokes' in the plot, and rows are individual entries, eg neurons

    subsets:    List
                A list of skids to subset the data by when plotting. If a list of lists is given, will overlay multiplpe plots on the axis, for each subset.

    size        Tuple
                Defaults to (10,10). With a larger number of axis, these figures tend to get very crouded, so big is good.

    Returns
    -------

    plot        fig, ax
                pyplot fig and axis objects

    """

    # Set up figure and Axis:
    # sort the figure size:
    fig = plt.figure(figsize = size)


    # Determine how many number of spokes we will have
    categories = list(data.columns)
    N = len(categories)
    # Work out at whhat angle eac spoke has to be, and make sure they form a full circle
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]#

    # Initialise the axis
    ax = plt.subplot(111, polar = True)
    # Draw spokes
    plt.xticks(angles[:-1],categories)

    # if only one data set:
    if subsets is None:
        # Add data
        # Take sum of each row over the total sum of all rows
        summed = list(data.sum()/data.sum().sum())
        summed += summed[:1]
        # Plot on axis
        ax.plot(angles, summed)
        ax.fill(angles, summed, alpha = 0.1)
    # if more than one subset:
    elif sum([isinstance(subsets[0],list) for i in subsets]) > 1:
        for subset in subsets:
            # Add data
            # Take sum of each row over the total sum of all rows
            summed = list(data.loc[subset].sum()/data.loc[subset].sum().sum())
            summed += summed[:1]
            # Plot on axis
            ax.plot(angles, summed)
            ax.fill(angles, summed, alpha = 0.1)
    # if only one subset...
    else:
        data = data.loc[subsets]
        # Add data
        # Take sum of each row over the total sum of all rows
        summed = list(data.sum()/data.sum().sum())
        summed += summed[:1]
        # Plot on axis
        ax.plot(angles, summed)
        ax.fill(angles, summed, alpha = 0.1)
