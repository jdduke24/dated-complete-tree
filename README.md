# dated-complete-tree
Take the Open Tree of Life, resolve its topology, and interpolate dates for all nodes.

Resolution can be performed repeatedly to generate a distribution of dated complete trees.

Optionally, compute a distribution of evolutionary distinctiveness scores for each leaf node.

## Prerequisites:

- The Open Tree of Life. By default we look for a folder ./opentree16.1_tree containing annotations.json and labelled_supertree/labelled_supertree_ottnames.tre
- The Open Tree Taxonomy. By default we look for a folder called ./ott3.7.3 containing taxonomy.tsv.
- The Python library chronosynth, available at https://github.com/OpenTreeOfLife/chronosynth/. This should be available in your PYTHONPATH.
- A folder to cache the chronosynth output. By default we use ./chronosynth_date_info.

## Usage:

There are three options:

1. The file main.py can be run from the command line. For available options run `python main.py --help`. For example:
> `python main.py --num_trees=10 --pd_clades=pd_clades.txt`
will produce 10 trees with different topologies, and output a text file with phylogenetic diversity distributions for the clades specified in your text file ./pd_clades.txt (one Open Tree node name per line).

2. The file main_non_exec.py contains code you can run from your favourite Python IDE.

3. A Jupyter notebook main_edge2.ipynb is included, which will open a Jupyter (IPython) notebook to step through the process of loading trees and computes EDGE2 scores.
