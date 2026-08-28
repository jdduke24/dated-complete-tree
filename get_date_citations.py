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


"""Get citations for the trees used in node_ages.json.
Expects a local clone of the phylesystem repository - adjust the path in line 70.
"""

import json

dates = json.load(open('date_cache/node_ages.json'))

my_dates = {}
my_sources = {}
for node in dates['node_ages']:
    my_dates[node] = []
    my_sources[node] = []
    for dat in dates['node_ages'][node]:
        my_dates[node].append(float(dat['age']))
        my_sources[node].append(dat['source_id'])

all_source_trees = set()
for node in my_sources:
    for source_tree in my_sources[node]:
        all_source_trees.add(source_tree)


seen_trees = set()
refs_text = []
for source_tree in all_source_trees:
    file_name = source_tree.split('@')[0]

    if file_name in seen_trees:
        continue
    else:
        seen_trees.add(file_name)

    phyle_folder = file_name[:3] + file_name[-2:]

    tree_name = source_tree.split('@')[1]

    print(source_tree, phyle_folder, file_name, tree_name)

    tree_dict = json.load(open("../phylesystem-1/study/%s/%s/%s.json" % (phyle_folder,
                                                                         file_name,
                                                                         file_name),
                                                                         'r'))

    ref_text = tree_dict["nexml"]["^ot:studyPublicationReference"].strip()
    ref_text = ref_text.replace('\n', '')

    refs_text.append(ref_text + '\n')


refs_text.sort()

with open('date_source_refs.txt', 'wt') as fout:
    fout.writelines(refs_text)
