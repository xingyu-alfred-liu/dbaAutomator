import os

bader_path = "/home/xingyu/papernow/PAHs-sf-silicon/PFA/homo-lumo/CARREU/bader"
HOMO_charge_path = "/home/xingyu/software/dba_automation/data/HOMO.cube"

def run_dba(bader_path, charge_path):
    os.system(bader_path+' '+charge_path)


if __name__ == "__main__":
    run_dba(bader_path, HOMO_charge_path)
