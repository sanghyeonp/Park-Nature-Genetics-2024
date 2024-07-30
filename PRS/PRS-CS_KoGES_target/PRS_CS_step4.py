import csv
import code
# code.interact(local=dict(globals(), **locals()))
import os

file = "score_sum_indiv.profile"
dir, _ = os.path.split(file)

with open(file, 'r') as f, open(file + '.reformat', 'wt', newline='') as f2, open(os.path.join(dir, "PRS_MetS_noUKB.KoGES.csv"), 'wt', newline='') as f3:
    reader = csv.reader(f, delimiter=' ')
    writer = csv.writer(f2, delimiter='\t')
    writer2 = csv.writer(f3, delimiter=',')
    header = reader.__next__()
    writer.writerow([l for l in header if l])
    writer2.writerow(['eid', 'PRS'])
    for row in reader:
        new_row = [d for d in row if d]
        writer.writerow(new_row)
        writer2.writerow([new_row[1], new_row[-1]])
        # code.interact(local=dict(globals(), **locals()))