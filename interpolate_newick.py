# BSD 3-Clause License

# Copyright (c) 2025, Jonathan David Duke

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import argparse
import ete4
import sys

from dated_complete_tree import tree_dating


def interpolate_dates(args):
    tre = ete4.Tree(args.phylo, parser=1)

    tree_dating.assign_dates_from_file(tre, args.ages)

    if args.use_bladj_consistency_fix:
        tree_dating.remove_inconsistent_dates(tre)
    else:
        tree_dating.dq_date_removal(tre)

    tree_dating.date_labelling(tre)

    if args.algo == "EQS-L":
        tree_dating.impute_missing_dates(tre)
    elif args.algo == "EQS-S":
        tree_dating.impute_missing_dates(tre, l=0)
    elif args.algo == "LnN":
        tree_dating.impute_missing_dates(tre, use_logN_model=True)
    elif args.algo == "BM":
        import numpy as np
        date_interpolation_rng = np.random.default_rng(seed=1)
        tree_dating.impute_missing_dates(tre, use_birth_model=True, rng=date_interpolation_rng)
    else:
        tree_dating.impute_missing_dates(tre, l=0.25)

    tree_dating.compute_branch_lengths(tre)

    if args.output_file:
        tree_dating.write_tree_with_branch_lengths(tre, args.output_file)
    else:
        from ete4.parser.newick import DIST, NAME
        MY_DIST = DIST.copy()
        MY_DIST['write'] = lambda x: "%.7f" % float(x)
        custom_parser = {
            'leaf':     [NAME, MY_DIST],
            'internal': [NAME, MY_DIST]
        }

        sys.stdout.write(tre.write(parser=custom_parser, format_root_node=True))


def main():
    parser = argparse.ArgumentParser(
        description=(
            """Interpolate missing dates (ages) on a given newick phylogenetic tree. Known ages should be listed
            in the given ages file. Leaf nodes without a date are assumed to have a date of 0. Must include an age on the root node."""
        )
    )

    parser.add_argument("--phylo",
                        help="Name of Newick-format text file specifying the tree topology. Defaults to 'phylo', like BLADJ.",
                        default="phylo")

    parser.add_argument("--ages",
                        help="Tab-separated files with two columns: node name and node ages. Defaults to 'ages', like BLADJ.",
                        default="ages")

    parser.add_argument("--output_file",
                        help="If specified, will write out a Newick file with branch lengths to this filename. Otherwise, will write to stdout, like BLADJ.",
                        default=None)

    parser.add_argument("--algo",
                        help="Interpolation algorithm to use. Defaults to EQS-LS, our recommended algorithm. To ape BLADJ behaviour, use EQS-L.",
                        choices=["EQS-L", "EQS-S", "EQS-LS", "LnN", "BM"],
                        default="EQS-LS")

    parser.add_argument("--use_bladj_consistency_fix",
                        help="""If this flag is used, the node ages will be made consistent using a top-down algorithm as in the BLADJ implementation. If absent, we use
                        our improved consistency algorithm that can keep more ages.""",
                        action="store_true")

    args = parser.parse_args()
    interpolate_dates(args)


if __name__ == "__main__":
    main()
