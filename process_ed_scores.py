import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

all_scores = []

fin = open("output/dated_tree_ed_scores.txt", "r")
for i, line in enumerate(fin):
    if i == 0:
        continue

    parts = line.split('\t')

    all_scores.append((float(parts[11]), float(parts[9]), float(parts[13]), parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]))
fin.close()

all_scores.sort(reverse=True)

fout = open("top_ed.csv", "w")
for i in range(len(all_scores)):
    if i > 100000:
        break

    fout.write("%s,%s,%s,%s,%s,%s,%s,%f,%f,%f\n" % (all_scores[i][3],
                                                    all_scores[i][4],
                                                    all_scores[i][5],
                                                    all_scores[i][6],
                                                    all_scores[i][7],
                                                    all_scores[i][8],
                                                    all_scores[i][9],
                                                    all_scores[i][0],
                                                    all_scores[i][1],
                                                    all_scores[i][2]))

fout.close()
