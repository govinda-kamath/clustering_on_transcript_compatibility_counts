import os
import multiprocessing as mp

hisat_idx_path = '/data/SS_RNA_seq/Zeisel/reference_transcriptome/Mus_musculus.GRCm38.rel79.cdna.all.HISAT'

def run_hisat(cell_tuple):
    cell,hisat_idx_path = cell_tuple[0],cell_tuple[1]
    outpath = './hisat/'+cell.split('/')[-1].replace('.fastq.gz','.sam')
    hisat_cmd = 'hisat --no-spliced-alignment -x '+hisat_idx_path+' -U '+cell+' -S '+outpath
    print hisat_cmd
    os.system(hisat_cmd)

os.system('mkdir -p ./hisat/')
reads_path = '/data/SS_RNA_seq/Zeisel/cleaned_reads_subsample100/'
cells = [(reads_path+cell,hisat_idx_path) for cell in os.listdir(reads_path) if 'fastq' in cell]

pool = mp.Pool(processes = 40)
pool.map(run_hisat,cells)
