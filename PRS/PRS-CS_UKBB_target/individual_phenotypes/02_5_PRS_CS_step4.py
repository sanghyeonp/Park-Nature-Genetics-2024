import csv
import os
import code

OUT="./PRS"
Traits=['BMI', 'WC', 'T2D', 'FG', 'HTN', 'TG', 'HDL']

for trait in Traits:
    file = os.path.join(OUT, trait, "score_sum_indiv.profile")
    
    with open(file, 'r') as f, open(file + '.reformat', 'wt', newline='') as f2:
        reader = csv.reader(f, delimiter=' ')
        writer = csv.writer(f2, delimiter='\t')
        for row in reader:
            new_row = [d for d in row if d]
            writer.writerow(new_row)