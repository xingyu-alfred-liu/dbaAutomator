import os

bader_path = "/home/xingyu/papernow/PAHs-sf-silicon/PFA/homo-lumo/CARREU/bader"
HOMO_charge_path = "/home/xingyu/software/dba_automation/data/"
HOMO_charge_file = HOMO_charge_path + 'HOMO.cube'
output_file = HOMO_charge_path + 'output'

def run_dba(bader_path, charge_file, output_file):
    os.system(bader_path+' '+charge_file+' > '+output_file)


if __name__ == "__main__":
    run_dba(bader_path, HOMO_charge_file, output_file)
    os.system('mv *dat ../data')
