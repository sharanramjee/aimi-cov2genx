import pickle, tempfile, subprocess
from threading import Thread

NETMHCPAN_LOCATION = "data/netMHCpan"
DICTS_FOLDER = "data/dicts"
MHCS = ["HLA-A01:01", "HLA-A02:01"]#, "HLA-A03:01", "HLA-A11:01", "HLA-A24:02", "HLA-B40:01", "HLA-C01:02", "HLA-C04:01", "HLA-C07:01", "HLA-C07:02"]
GENES = ['e']#, 'm', 'n', 's']
MIN_PEPTIDE_LEN = 8
MAX_PEPTIDE_LEN = 10
NETMHCPAN_COLUMNS = ["pos","mhc","peptide","core","of","gp","gl","ip","il","icore","identity","score_el","bindlevel"]
DICT_ENTRY_SIZE = 4
def parse_netmhcpan_output(netmhcpan_output: str, expected_count: int) -> [{str: str}]:
    netmhcpan_output = netmhcpan_output.split('\n')[49:]
    result = []
    for line in netmhcpan_output:
        line = [entry for entry in line.split() if entry != ""][:len(NETMHCPAN_COLUMNS)]
        if len(line) != 13:
            break
        result.append({key: value for key, value in zip(NETMHCPAN_COLUMNS, line)})
    if len(result) != expected_count:
        raise ValueError("Invalid number of peptides parsed from netmhcpan")
    return result
print(','.join(NETMHCPAN_COLUMNS + [f"dicts_column_{i}" for i in range(DICT_ENTRY_SIZE)] + ['gene']))
for gene in GENES:
    with open(f"{DICTS_FOLDER}/{gene}_dict.pickle", 'rb') as pickled_dict:
        dictionary = pickle.load(pickled_dict)
    netmhcpan_output = {mhc: "" for mhc in MHCS}
    peptides_count = 0
    with tempfile.NamedTemporaryFile('w') as netmhcpan_input:
        for peptide in dictionary:
            if len(peptide) >= MIN_PEPTIDE_LEN and len(peptide) <= MIN_PEPTIDE_LEN:
                peptides_count += 1
                print(peptide, file=netmhcpan_input)
        netmhcpan_input.flush()
        def run_netmhcpan(mhc: str) -> None:
            netmhcpan_output[mhc] = parse_netmhcpan_output(subprocess.check_output(args=[NETMHCPAN_LOCATION, "-a", mhc, "-p", "-f", netmhcpan_input.name]).decode("utf-8"), peptides_count)
        threads = [Thread(target=run_netmhcpan, args=(mhc,)) for mhc in MHCS]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
    for mhc, peptide_entries in netmhcpan_output.items():
        for peptide_entry in peptide_entries:
        	for dictionary_entry in dictionary[peptide_entry["peptide"]]:
                    for column in NETMHCPAN_COLUMNS:
                        print(peptide_entry[column], end=',')
                    for i in range(DICT_ENTRY_SIZE):
                        print(dictionary_entry[i], end=',')
                    print(gene.upper() + "_gene")
