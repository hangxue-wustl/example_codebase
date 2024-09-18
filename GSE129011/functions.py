import os
import subprocess

def trinity_helper(input_folder, library, seqtype):
    outdir = f"{input_folder}/trinity_{library}_outdir"
    os.makedirs (outdir, exist_ok=True)
    res_out = f"{outdir}/{library}_trinity.res"
    std_err = f"{outdir}/{library}_trinity.log"
    res_out_file = open(res_out, 'w')
    std_err_file = open(std_err, 'w')

    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "--bind", f"{input_folder}/trinity_{library}_outdir:/outdir",
            "/global/scratch/users/hangxue/sif/trinityrnaseq_latest.sif",
            "Trinity", 
            "--seqType", f"{seqtype}",
            "--single",f"/input/{library}.gz",
            "--output", "/outdir/trinity",
            "--CPU","16",
            "--min_contig_length", "100",
            "--max_memory", "20G"]

    subprocess.run(cmd,stdout=res_out_file, stderr = std_err_file)
    res_out_file.close()
    std_err_file.close()
    # trinitystat(library, outdir)
    print (outdir)

def trinitystat(library, outdir):
    res_out = f"{outdir}/{library}_trinitystats.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{outdir}:/outdir",
            "/global/scratch/users/hangxue/sif/trinityrnaseq_latest.sif",
            "/usr/local/bin/util/TrinityStats.pl",
            "/outdir/Trinity.fasta"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def minimap2_index(input_folder, fasta):
    res_out = f"{input_folder}/{fasta}_minimap2.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "/global/scratch/users/hangxue/sif/minimap2_latest.sif", 
            "minimap2", 
            "-d", f"/input/{fasta}.mmi",
            f"/input/{fasta}"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def minimap2_map(input_folder1, input_folder2, reference, trinity_assembly):
    outdir = f"{input_folder2}/minimap2_outdir"
    os.makedirs (outdir, exist_ok=True)
    res_out = f"{outdir}/minimap2.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder1}:/input1",
            "--bind", f"{input_folder2}:/input2",
            "--bind", f"{outdir}:/outdir",
            "/global/scratch/users/hangxue/sif/minimap2_latest.sif", 
            "minimap2", 
            "-a",f"/input1/{reference}",f"/input2/{trinity_assembly}",
            "-o", f"/outdir/{trinity_assembly}_minimap.sam"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def blast_makedb(input_folder, reference):
    res_out = f"{input_folder}/{reference}_blastmakedb.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "/global/scratch/users/hangxue/sif/blast_latest.sif", 
            "makeblastdb", 
            "-in", f"/input/{reference}",
            "-out", f"{reference}_db", 
            "-parse_seqids",
            "-dbtype", "nucl"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def blastn(input_folder, query, db, outdir):
    res_out = f"{input_folder}/{reference}_blastmakedb.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "--bind", f"{outdir}:/output",
            "/global/scratch/users/hangxue/sif/blast_latest.sif", 
            "blastn", 
            "-query", f"/output/{query}",
            "-db", f"/input/{db}", 
            "-outfmt", "6 qseqid qstart qend sstart send pident length mismatch gapopen evalue bitscore",
            "-out", f"/output/trinity_blast.out"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def samtools_sort (input_folder, file):
    res_out = f"{input_folder}/{file}.sort"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "/global/scratch/users/hangxue/sif/samtools_latest.sif", 
            "samtools", 
            "sort", 
            f"/input/{file}"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def samtools_cov (input_folder, file):
    res_out = f"{input_folder}/{file}_coverage.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "/global/scratch/users/hangxue/sif/samtools_latest.sif", 
            "samtools", 
            "coverage", 
            f"/input/{file}"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)

def samtools_dep(input_folder, file):
    res_out = f"{input_folder}/{file}_dep.res"
    res_out_file = open(res_out, 'w')
    cmd = ["singularity" ,"exec" , 
            "--bind", f"{input_folder}:/input",
            "/global/scratch/users/hangxue/sif/samtools_latest.sif", 
            "samtools", 
            "depth", 
            f"/input/{file}"]
    subprocess.run(cmd, stdout = res_out_file, stderr = True)