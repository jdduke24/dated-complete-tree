# dated-complete-tree
Code that takes the Open Tree of Life, resolves all polytomies, and interpolates dates for all nodes. For a full explanation of the algorithms and the results of using them, see [this preprint](https://doi.org/10.64898/2026.03.05.709771 ).

Resolution and interpolation can be performed repeatedly to generate distributions of dated complete trees. Dates can be sampled from available date sources in the [Open Tree phylesystem](github.com/OpenTreeOfLife/phylesystem-1).

Optionally, you can also compute a distribution of evolutionary distinctiveness scores (you can choose ED/EDGE or ED2/EDGE2) for each leaf node.

Pre-computed median trees and tree distributions can be found [at the accompanying Zenodo dataset](https://doi.org/10.5281/zenodo.19049120).

## Prerequisites:

- The Open Tree of Life. By default we look for a folder `./opentree16.1_tree/` containing the files `annotations.json` and `labelled_supertree/labelled_supertree_ottnames.tre`.

- The Open Tree Taxonomy. By default we look for a folder called `./ott3.7.3/` containing `taxonomy.tsv`.

- The Python library chronosynth, available at https://github.com/OpenTreeOfLife/chronosynth/. This should be available in your PYTHONPATH.

- A folder to cache the chronosynth output. By default we use `./chronosynth_date_info/`.

## Usage:

There are three options:

1. The file main.py can be run from the command line. For available options run:
`python main.py --help`.
For example:
 `python main.py --num_trees=10 --pd_clades=pd_clades.txt`
will produce 10 trees with different topologies and a text file with phylogenetic diversity (PD) distributions for the clades specified in your text file `./pd_clades.txt` (one Open Tree node name per line, e.g. Eukaryota_ott304358). The default output folder for the Newick-format trees and the PD file is `./output/`.

2. The file `main_non_exec.py` contains code you can edit and run from your favourite Python IDE.

3. A Jupyter notebook `edge2_notebook.ipynb` is included, which will open a Jupyter (IPython) notebook to step through the process of loading trees and computes EDGE2 scores.

## Subtrees

If you want to work with only a subtree or subset of species, see the notebook `getting_a_subtree.ipynb`.


## Reproducing the pre-computed tree distributions

Should you wish to reproduce the distributions of trees used in the paper, you can use the following commands:

- equal_splits_topo.tar.gz: `python main.py --num_trees=501 --pd_clades=pd_clades.txt`

- equal_splits_both.tar.gz: `python main.py --num_trees=501 --num_date_samples=3 --pd_clades=pd_clades.txt`

- birth_model_topo.tar.gz: `python main.py --use_birth_model --num_trees=101 --pd_clades=pd_clades.txt`

- birth_model_both.tar.gz: `python main.py --use_birth_model --num_trees=101 --num_date_samples=3 --pd_clades=pd_clades.txt`