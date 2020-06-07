import csv
import pandas as pd


def modify_row(row):
    try:
        peptide = list(row[8])
        if len(peptide) == 8:
            peptide += ['-', '-']
        elif len(peptide) == 9:
            peptide += ['-']
        elif len(peptide) > 10:
            peptide = peptide[:10]
        pos = int(row[2])
        edit_dist = int(row[17])
        hla = row[0]
        location = row[1].split('/')[1]
        date = row[1].split('|')[-1]
        date_int = int(''.join(date.split('-')))
        l_score = float(row[3])
        r_score = float(row[4])
        mid_score = 1-float(row[5])
        existence_score = l_score * mid_score * r_score
        binding_score = float(row[16])
        ma_binding_score = existence_score * binding_score
        new_row = peptide + [pos] + [edit_dist] + [hla] + [location] + [date_int] + [ma_binding_score]
        return new_row
    except:
        print(row)


def data_generator(filename):
    try:
        with open(filename, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            next(csv_reader)
            while True:
                yield(next(csv_reader))
                for j in range(10000):
                    next(csv_reader)
    except:
        return


def make_dataframe():
    new_fields = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'pos', 'edit_dist', 'hla', 'location',
                  'days_since_g0', 'ma_binding_score']
    df = pd.DataFrame(columns=new_fields)
    generator = data_generator('peptides_all.csv')
    row_idx = 0
    for row in generator:
        df.loc[row_idx] = modify_row(row)
        row_idx += 1
        print(row_idx)
    return df
