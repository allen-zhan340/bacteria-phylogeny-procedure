import os
import sys
import datetime
import argparse
import time

today = datetime.date.today().strftime('Results_%b%d')


def get_parameters():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', required=True, action='store', help='Input a directory contains genome files, and the genomes must end with (.fa, .fna., .fasta)')
    parse.add_argument('-a', required=True, action='store', help='Output directory contains proteins (.faa)')
    parse.add_argument('-t', required=True, action='store', help='Threads used for prokka,orthofinder')
    args = parse.parse_args()
    return args



def prokka_(genomes_dir, threads):
    endslist = ['fa', 'fna', 'fasta']
    to_dir = 'final_prokka'
    if not os.path.exists(to_dir):
        os.mkdir(to_dir)
    for i in os.listdir(genomes_dir):
        ends_ = i.split(".")[-1]
        if ends_ in endslist:
            genomes_path = os.path.join(genomes_dir, i)
            outdir = f'{to_dir}/{i.split(".f")[0]}'
            prefix = f'{i.split(".f")[0]}'
            cmd = f"prokka --outdir {outdir} --prefix {prefix} --noanno --norrna --notrna --locustag '{prefix}|ORF' --cpus {threads} {genomes_path}"
            print(cmd)
            os.system(cmd)

def orthofinder_(faa_dir, threads):
    if not os.path.exists(faa_dir):
        os.mkdir(faa_dir)
    try:
        cp_cmd = f'cp final_prokka/*/*faa {faa_dir}'
        os.system(cp_cmd)
    except Exception:
        pass
    if not os.path.exists(f'{faa_dir}/OrthoFinder'):
        cmd_ = f'orthofinder -og -f {faa_dir} -t {threads}'
        os.system(cmd_)


def clusto_(faa_dir):
    out_put = 'final_clustalo_fasta'
    if not os.path.exists(out_put):
        os.mkdir(out_put)
    #for i in os.listdir("faa_dir/OrthoFinder/"):
	
    for i in os.listdir("faa_dir/OrthoFinder/" + today + "/Single_Copy_Orthologue_Sequences"):
        name = i.split(".f")[0]
        to_name = f'{out_put}/{name}.fasta'
        cmd = f'clustalo -i {faa_dir}/OrthoFinder/{today}/Single_Copy_Orthologue_Sequences/{i} -o {to_name}'
        os.system(cmd)
        print(cmd)

def align_():
    cmd = 'python AlignConcat.py -i final_clustalo_fasta -o final.aln'
    os.system(cmd)


def gblock_():
    cmd = 'sh Gblocks.sh final.aln'
    os.system(cmd)


def iqtree_():
    cmd = 'iqtree -s final.aln-gb -m LG+F+R4 -bb 1000 -nt 45'
    os.system(cmd)

if __name__ == '__main__':
    #python tree_build_orthofinder.py -i genomes -a faa_dir -t 60 

    arge = get_parameters()
    genomes_dir = arge.i
    faa_dir = arge.a
    threads = arge.t


    clusto_dir = 'final_clustalo_fasta'

    if not os.path.exists(faa_dir):
        prokka_(genomes_dir, threads)

    
    orthofinder_(faa_dir, threads)
    time.sleep(5)
 
    if not os.path.exists(clusto_dir):
        clusto_(faa_dir)
    if not os.path.exists('final.aln'):
        align_()

    #if not os.path.exists('final.aln-gb'):
    gblock_()

    iqtree_()






