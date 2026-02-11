"""
Phylogenetic tree plotter using matplotlib and ete4
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from ete4 import Tree


def name_to_simple_name(name):
    parts = name.split("_")
    ret_str = parts[0]
    if len(parts) > 0:
        for p in parts[1:-1]:
            ret_str = ret_str + " " + p
    return ret_str.strip()


def digit_to_hex(dig):
    """Convert an integer between 0 and 15 to a string representing a hex value between 0 and F."""
    if dig < 10:
        return str(dig)
    elif dig == 10:
        return "A"
    elif dig == 11:
        return "B"
    elif dig == 12:
        return "C"
    elif dig == 13:
        return "D"
    elif dig == 14:
        return "E"
    elif dig == 15:
        return "F"
    else:
        print(dig)
        raise ValueError


def pct_to_color_hex_str(pct):
    """Convert a number between 0 and 1 to a string representing a hex value between 00 and FF."""
    p = pct * 255
    first = int(p/16)
    second = int(p % 16)
    return digit_to_hex(first) + digit_to_hex(second)


def get_tree_layout(tree):
    """
    Calculate x,y coordinates for all nodes in the tree.
    Returns a dictionary mapping nodes to (x, y) coordinates.
    """
    positions = {}

    # First pass: assign y-coordinates to leaves (sequential numbering)
    leaf_y_counter = [0]  # Use list to allow modification in nested function

    def assign_positions(node, depth):
        """
        Recursively calculate positions for all nodes.
        depth: cumulative distance from root
        """
        if node.is_leaf:
            # Leaves get sequential y-positions
            y = leaf_y_counter[0]
            leaf_y_counter[0] += 1
            x = depth - (node.dist if node.dist else 0)
            positions[node] = (x, y)
            return y
        else:
            # Process children first
            current_depth = depth - (node.dist if node.dist else 0)
            child_y_positions = []

            for child in node.children:
                child_y = assign_positions(child, current_depth)
                child_y_positions.append(child_y)

            # Internal node y-position is the average of its children
            y = sum(child_y_positions) / len(child_y_positions)
            x = depth - (node.dist if node.dist else 0)
            positions[node] = (x, y)
            return y

    # Start from root with depth 0
    assign_positions(tree, 4267.7+tree.dist)

    return positions


def plot_tree(tree, cmap, ax=None, figsize=(10, 8), align_labels=True):
    """
    Plot a phylogenetic tree using matplotlib.

    Parameters:
    -----------
    tree : ete4.Tree
        The tree object to plot
    ax : matplotlib axes, optional
        Axes to plot on. If None, creates new figure
    figsize : tuple
        Figure size if creating new figure

    Returns:
    --------
    fig, ax : matplotlib figure and axes
    """

    to_name = set(["cellular",
                   "Eukaryota",
                   "Opisthokonta",
                   "Fungi",
                   "Holozoa",
                   "Metazoa",
                   "Bilateria",
                   "SAR",
                   "Chloroplastida",
                   "Ecdysozoa"])

    def label_name(node_name):
        for name in to_name:
            if name in node_name:
                return True

        return False


    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Calculate positions for all nodes
    positions = get_tree_layout(tree)

    all_pds = []
    all_x = []
    for node in tree.traverse(strategy="preorder"):
        all_pds.append(node.props["pd"])
        x_node, y_node = positions[node]
        all_x.append(x_node)

    max_pd = np.log(np.max(all_pds))
    min_pd = np.log(np.min(all_pds))

    min_x = np.log10(np.min(all_x))

    # Draw edges
    for idx, node in enumerate(tree.traverse(strategy="levelorder")):
        # if node.is_root:
        #     pd_conversion = (np.log(node.props["pd"]) - min_pd) / (max_pd - min_pd)
        #     line_width = 2 + pd_conversion * 5
        #     node.add_prop("line_width", line_width)
        #     continue

        # Get parent and current node positions
        x_node, y_node = positions[node]
        if node.up:
            parent = node.up
            x_parent, y_parent = positions[parent]
        else:
            x_parent = x_node + node.dist
            y_parent = y_node

        x_parent = np.log10(x_parent)
        x_node = np.log10(x_node)

        # Get line properties from props dictionary if available
        props = getattr(node, "props", {})
        # line_color = props.get("color", "black")

        pd_conversion = (np.log(node.props["pd"]) - min_pd) / (max_pd - min_pd)

        line_width = 1 + pd_conversion*0.8
        # line_width = 2
        line_style = props.get("linestyle", "-")

        node.add_prop("line_width", line_width)

        rgba = cmap(pd_conversion)
        color_str = "#"
        for i in range(3):
            color_str += pct_to_color_hex_str(rgba[i])
        line_color = color_str

        # Draw horizontal line
        z = 5*idx
        node.add_prop("zorder", z)
        ax.plot([x_parent, x_node], [y_node, y_node],
                color=line_color, linewidth=line_width, linestyle=line_style, zorder=z)

        # print(node.name, z)

        # Draw vertical line
        if node.up:
            z = node.up.props["zorder"]
            if y_node > y_parent:
                # print("  down:", z-7)
                ax.plot([x_parent, x_parent], [y_parent, y_node],
                        color=line_color, linewidth=line_width, linestyle=line_style, zorder=z-7)
                ax.plot([x_parent, x_parent], [y_parent, y_node],
                        color="white", linewidth=line_width*1.3, linestyle=line_style, zorder=z-8)
            else:
                # print("  up:", z-2)
                ax.plot([x_parent, x_parent], [y_parent, y_node],
                        color=line_color, linewidth=line_width, linestyle=line_style, zorder=z-1)
                ax.plot([x_parent, x_parent], [y_parent-0.1, y_node],
                        color="white", linewidth=line_width*1.3, linestyle=line_style, zorder=z-2)

    # Draw leaf labels
    if not align_labels:
        for leaf in tree.leaves():
            x, y = positions[leaf]
            x = np.log10(x)

            ax.text(x*0.999, y, name_to_simple_name(leaf.name),
                    verticalalignment="center",
                    horizontalalignment="right",
                    fontsize="small",
                    fontfamily="Arial",
                    color="black")
    else:
        # pd_conversion = (np.log(tree.props["pd"]) - min_pd) / (max_pd - min_pd)
        # rgba = cmap(pd_conversion)
        # color_str = "#"
        # for i in range(3):
        #     color_str += pct_to_color_hex_str(rgba[i])

        # x, y = positions[tree]
        # x = np.log10(x)
        # # square at rootnode
        # squ_pd_conv = 0.45*pd_conversion + 0.55

        # rect_x_size = 0.008
        # rect_y_size = rect_x_size*50

        # rect = patches.Rectangle(((x-rect_x_size)+(rect_x_size-rect_x_size*squ_pd_conv/2),
        #                           (y-rect_y_size/2)+(rect_y_size-rect_y_size*squ_pd_conv)/2),
        #                           rect_x_size*squ_pd_conv,
        #                           rect_y_size*squ_pd_conv,
        #                           linewidth=1,
        #                           edgecolor=color_str,
        #                           facecolor=color_str,
        #                           zorder=10000)

        # ax.add_patch(rect)

        for leaf in tree.leaves():
            x, y = positions[leaf]
            x = np.log10(x)

            ax.text(min_x-0.015, y, name_to_simple_name(leaf.name),
                    verticalalignment="center",
                    horizontalalignment="right",
                    fontsize="x-small",
                    fontfamily="Arial",
                    color="black")

            ax.plot([min_x-0.01, x], [y, y], color="#cfcfcf", linewidth=0.6, linestyle="dashed", dashes=(6,3), zorder=-1000)

            pd_conversion = (np.log(leaf.props["pd"]) - min_pd) / (max_pd - min_pd)
            rgba = cmap(pd_conversion)
            color_str = "#"
            for i in range(3):
                color_str += pct_to_color_hex_str(rgba[i])

            pd_conversion = 0.45*pd_conversion + 0.55

            rect_x_size = 0.02
            rect_y_size = 1

            # square at node
            # rect = patches.Rectangle(((min_x-0.026)+(rect_x_size-rect_x_size*pd_conversion/2),
            #                           (y-rect_y_size/2)+(rect_y_size-rect_y_size*pd_conversion)/2),
            #                           rect_x_size*pd_conversion,
            #                           rect_y_size*pd_conversion,
            #                           linewidth=0,
            #                           edgecolor=color_str,
            #                           facecolor=color_str)

            # ax.add_patch(rect)

            # triangle at left
            # from matplotlib.path import Path

            # x_tri = (min_x-0.026)+(rect_x_size-rect_x_size*pd_conversion/2)
            # y_tri = (y-rect_y_size/2)+(rect_y_size-rect_y_size*pd_conversion)/2

            # path = Path([[x_tri, y_tri],
            #              [x_tri, y_tri+rect_y_size*pd_conversion],
            #              [x_tri+rect_x_size*pd_conversion, y_tri+rect_y_size*pd_conversion/2],
            #              [x_tri,y_tri]])

            # tri = patches.PathPatch(path, linewidth=0, color=color_str)

            # ax.add_patch(tri)

            # triangle at node
            from matplotlib.path import Path

            x_tri = (x-rect_x_size)+(rect_x_size-rect_x_size*pd_conversion/2)
            y_tri = (y-rect_y_size/2)+(rect_y_size-rect_y_size*pd_conversion)/2

            path = Path([[x_tri, y_tri],
                         [x_tri, y_tri+rect_y_size*pd_conversion],
                         [x_tri+rect_x_size*pd_conversion, y_tri+rect_y_size*pd_conversion/2],
                         [x_tri,y_tri]])

            tri = patches.PathPatch(path, linewidth=0, color=color_str)

            ax.add_patch(tri)

            # repeat square at node
            # rect = patches.Rectangle(((x-0.014)+(rect_x_size-rect_x_size*pd_conversion/2),
            #                           (y-0.4)+(rect_y_size-rect_y_size*pd_conversion)/2),
            #                           rect_x_size*pd_conversion,
            #                           rect_y_size*pd_conversion,
            #                           linewidth=1,
            #                           edgecolor=color_str,
            #                           facecolor=color_str)

            # ax.add_patch(rect)

        # for node in tree.traverse():
        #     if node.props["imputed_date"]:
        #         continue

        #     x, y = positions[node]
        #     x = np.log10(x)

        #     # circle calibration nodes
        #     rect_x_size = 0.02
        #     rect_y_size = rect_x_size*50

        #     rect = patches.Ellipse((x,y),
        #                             rect_x_size*pd_conversion,
        #                             rect_y_size*pd_conversion,
        #                             fill=False,
        #                             rasterized=True,
        #                             linewidth=0.6,
        #                             edgecolor="black",
        #                             facecolor=None,
        #                             zorder=100000)

        #     ax.add_patch(rect)

    # draw important internal node labels
    for node in tree.traverse():
        if not node.is_leaf and label_name(node.name):
            x, y = positions[node]
            x = np.log10(x)

            tx = ax.text(x+0.0035, y-0.23, name_to_simple_name(node.name) if node is not tree else "All Life",
                   verticalalignment="top",
                   horizontalalignment="left",
                   fontsize="xx-small" if node is not tree else "x-small",
                   fontfamily="Arial",
                   color="black")

            # if "Bilateria" in node.name:
            #     tx.set_bbox(dict(facecolor="white", alpha=0.95, linewidth=0, boxstyle="round,pad=0.1"))

    # Clean up axes
    ax.set_xlabel("Time (Mya, log scale)", fontfamily="Arial", fontsize="x-small")
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_ylim((-1,72))

    yrs = list(range(250,4300,50))

    to_delete = []
    for i, x in enumerate(yrs):
        if (x > 500 and x % 100 != 0) or (x > 1000 and x % 200 != 0) or (x > 2000 and x % 400 != 0):
            to_delete.append(i)

    for i in to_delete[::-1]:
        del yrs[i]

    ax.set_xticks([np.log10(x) for x in yrs])
    ax.set_xticklabels([str(x) for x in yrs], fontfamily="Arial", fontsize="x-small")


    return fig, ax



# Create the example tree
# t = Tree("((A,B)internal, C)root;", parser=1)

# # Add props to some nodes for custom formatting
# for node in t.traverse():
#     if node.is_leaf:
#         if node.name == "A":
#             node.props = {"color": "red", "linewidth": 2.0, "label_color": "red"}
#         elif node.name == "B":
#             node.props = {"color": "blue", "linewidth": 2.0, "label_color": "blue"}
#         elif node.name == "C":
#             node.props = {"color": "green", "linewidth": 2.0, "label_color": "green"}
#     else:
#         # Color internal branches differently
#         if node.name == "internal":
#             node.props = {"color": "purple", "linewidth": 2.5}
#         elif node.name == "root":
#             node.props = {"color": "black", "linewidth": 1.5}
#         else:
#             node.props = {"color": "gray", "linewidth": 1.0}
#     node.dist = 1


import pickle
with open("figures/figure3/tree_test.pickle", "rb") as f:
    whole_tre = pickle.load(f)

whole_tre.dist = 30

for node in whole_tre.traverse(strategy="postorder"):
    if not node.is_leaf:
        if node.children[1].props["pd"] > node.children[0].props["pd"]:
            node.add_child(node.children[0].detach())


orig_cmap = plt.colormaps["plasma_r"]

import matplotlib.colors as colors
minval = 0.07
maxval = 1
n = 1000
cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=orig_cmap.name, a=minval, b=maxval),
        orig_cmap(np.linspace(minval, maxval, n)))

# Plot the tree
fig, ax = plot_tree(whole_tre, cmap, figsize=(7.5, 7.75), align_labels=True)



import matplotlib as mpl
c_map_ax = fig.add_axes([0.68, 0.09, 0.028, 0.35])
c_map_ax.axes.get_xaxis().set_visible(False)
# c_map_ax.axes.get_yaxis().set_visible(False)
cb = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cmap, orientation='vertical')

cb.set_label('Phylogenetic diversity (Byr, log scale)', fontfamily="Arial", fontsize="x-small")


# c_map_ax.set_label('Colour scale for PD per clade (Byr, log scale)')

all_pds = {}
pd_list = []
for node in whole_tre.traverse(strategy="preorder"):
    all_pds[name_to_simple_name(node.name)] = node.props["pd"]
    pd_list.append(node.props["pd"])
    print(node.name, node.props["date"], node.props["imputed_date"], node.props["pd"])

max_pd = np.log10(np.max(pd_list))
min_pd = np.log10(np.min(pd_list))

pds = 500*16**np.arange(0,5)

c_map_ax.yaxis.tick_left()
c_map_ax.yaxis.set_label_position("left")
c_map_ax.set_yticks((np.log10(pds) - min_pd) / (max_pd - min_pd))
c_map_ax.set_yticklabels([(str(round(x/1000,1)) if x < 1000 else str(int(x/1000))) for x in pds], fontfamily="Arial", fontsize="x-small")

named_pds = [("cellular organisms", "All Life"),
             ("Eukaryota", "Eukaryota"),
             ("Metazoa", "Metazoa"),
             ("Arthropoda", "Arthropoda"),
             ("Bacteria", "Bacteria"),
             ("Fungi", "Fungi"),
             ("Chloroplastida", "Chloroplastida"),
             ("Mollusca", "Mollusca"),
             ("Chordata", "Chordata"),
             ("Cnidaria", "Cnidaria"),
             ("Archaea", "Archaea"),
             ("Rhodophyta", "Rhodophyta"),
             ("Ambulacraria", "Ambulacraria"),
             ("Bryozoa", "Bryozoa"),
             ("Tardigrada", "Tardigrada"),
             ("Haptophyta", "Haptophyta"),
             ("Nematomorpha","Nematomorpha"),
             ("Gnathostomulida","Gnathostomulida"),
             ("Glaucophyta","Glaucophyta"),
             ("Placozoa", "Placozoa")
             ]

count = 0
for tree_name, my_name in named_pds:
    pd = all_pds[tree_name]
    y = (np.log10(pd) - min_pd) / (max_pd - min_pd)

    if count % 2 == 1:
        c_map_ax.plot([1, 1.5], [y, y], clip_on=False, linewidth=0.5, color="gray")
        c_map_ax.text(1.55, y+0.0001, "%s: %d Byr" % (my_name, round(pd/1000,0)), horizontalalignment="left", verticalalignment="center", fontfamily="Arial", fontsize="xx-small")
    else:
        c_map_ax.plot([1, 6.35], [y, y], clip_on=False, linewidth=0.5, color="gray")
        c_map_ax.text(6.4, y+0.0001, "%s: %d Byr" % (my_name, round(pd/1000,0)), horizontalalignment="left", verticalalignment="center", fontfamily="Arial", fontsize="xx-small")

    count += 1

fig.tight_layout()
fig.savefig("figures/figure3/mpl_tri_nodes7_narrow.tif", dpi=300)
fig.savefig("figures/figure3/mpl_tri_nodes7_narrow.pdf", dpi=300)
