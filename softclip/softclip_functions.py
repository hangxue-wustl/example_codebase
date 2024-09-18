import os
import glob
import shutil
import pysam
import pandas as pd
import numpy as np
import subprocess
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import time 
import multiprocessing as mp
from scipy.signal import savgol_filter, find_peaks, peak_widths
from statistics import mean 
            
def get_m_df():
    m_df = pd.read_csv('m_df.csv')
    m_df.index = m_df['Libraries']
    m_df.drop(columns = ['Libraries'], inplace = True)
    return (m_df)

def organize_folders(m_df):
    for library in m_df.index:
        os.makedirs(library, exist_ok = True)
        params = set_parameters(m_df,library)
        R1 = params.get("R1")
        R2 = params.get("R2")
        fasta = params.get("fasta")
        # in case file already exists in the correct folder
        if len(glob.glob(f"./{library}/{R1}")) < 1:
            R1 = "data/" + params.get("R1")
            shutil.copy(src = os.path.realpath(R1), dst = "./" + library)
        if len(glob.glob(f"./{library}/{R2}")) < 1:
            R2 = "data/" + params.get("R2")
            shutil.copy(src = os.path.realpath(R2), dst = "./" + library)
        if len(glob.glob(f"./{library}/{fasta}")) < 1:
            fasta =  "fasta/" + params.get("fasta")
            shutil.copy(src = os.path.realpath(fasta), dst = "./" + library)

def trimmomatic(m_df, clip_string):
    for library in m_df.index:
        params = set_parameters(m_df, library)
        R1 = library + "/" + params.get("R1")
        R2 = library + "/" +  params.get("R2")
        trimmed_R1 = params.get("trimmed_R1")
        trimmed_R2 = params.get("trimmed_R2")
        unpaired_R1 = params.get("unpaired_R1")
        unpaired_R2 = params.get("unpaired_R2")
        log = R1.replace(".fastq.gz","trim_out.log")
        working_dir = os.getcwd()
        if len(glob.glob(f"{trimmed_R1}")) < 1:
            com = ["docker","run","--rm",
                "-v", f"{working_dir}:/data",
                "--workdir", "/data",
                "trimmomatic:1.0.0",
                "trimmomatic","PE",
                R1, R2, 
                "-phred33",
                "-trimlog", log,
                trimmed_R1, unpaired_R1,
                trimmed_R2, unpaired_R2,
                clip_string,
                "SLIDINGWINDOW:4:30",
                "MINLEN:50"]
            logfile = open(f'{library}/{library}_trimmomatic.log', 'a')
            proc = subprocess.run(com, stdout=logfile, stderr=logfile)
            logfile.close()
            print("Trimmomatic has completed for " + library)
        else:
            print("Trimmomatic has previously completed for " + library)
    print("Trimmomatic has completed for all libraries")

def star_w_chimeric(m_df, n_proc):
    for library in m_df.index:
        star_w_chimeric_helper(m_df, library, n_proc)
    
def star_w_chimeric_helper(m_df, library, n_proc):
    working_dir = os.getcwd()
    params = set_parameters(m_df, library)
    reference_dir = params.get("reference_dir")
    reference_fa = params.get("reference_fa")
    star_out_prefix = params.get("star_out_prefix") + "CHIMERIC_"
    trimmed_R1 = params.get("trimmed_R1")
    trimmed_R2 = params.get("trimmed_R2")
    os.makedirs(reference_dir, exist_ok = True)
    shutil.copy(library + "/" + reference_fa,reference_dir)

    com = ["docker","run","--rm",
        "-v", f"{working_dir}:/data",
        "--workdir", "/data",
        "star:1.0.0",
        "STAR","--runThreadN", str(n_proc),
        "--genomeDir", reference_dir,
        "--chimOutType", "Junctions", "WithinBAM", "SoftClip",
        "--chimSegmentMin", "12", #here
        "--chimScoreJunctionNonGTAG", "0",
        "--limitBAMsortRAM 15000000000",
        "--outSAMtype BAM SortedByCoordinate",
        "--outSAMattributes NH HI AS nM MD",
        "--outFilterMismatchNoverLmax 0.3", #here
        "--alignIntronMax 1",
        #start here
        "--outSAMunmapped Within",
        "--twopassMode Basic",
        "--alignMatesGapMax 10000",
        "--chimJunctionOverhangMin 8",
        "--peOverlapNbasesMin 10",
        "--peOverlapMMp 0.1",
        "--chimMultimapNmax 20",
        #end here
        "--outFileNamePrefix", star_out_prefix,
        "--readFilesIn", trimmed_R1, trimmed_R2]
    subprocess.run(com,stdout= True, stderr = True)
    print("Star chimeric has completed for " + library)
    
def get_bam_list(m_df):
    bam_list = []
    for library in m_df.index:
        bam_list.append (glob.glob(f"{library}/{library}*CHIMERIC_Aligned.sortedByCoord.out.bam")[0])
    return (bam_list)

def samtools_index(file_list, n_cores):
    with mp.Pool(processes=n_cores) as p:
        p.map(samtools_index_helper, file_list)
    print("samtools_index has completed for all libraries")

def samtools_index_helper(bamname):
    working_dir = os.getcwd()
    cmd = ["docker" ,"run" , "--rm" ,
               "-v", f"{working_dir}:/data",
               "--workdir", "/data",
               "samtools:1.0.0",
               "samtools" ,"index", bamname]
    subprocess.run(cmd)
    return bamname

def samtools_flagstat(file_list, n_cores):
    with mp.Pool(processes=n_cores) as p:
        p.map(samtools_flagstat_helper, file_list)
    print("samtools_flagstat has completed for all libraries")

def samtools_flagstat_helper(bamname):
    working_dir = os.getcwd()
    cmd = ["docker" ,"run" , "--rm" ,
            "-v", f"{working_dir}:/data",
            "--workdir" , "/data",
            "samtools:1.0.0",
            "samtools" ,"flagstat",bamname]
    with open(f"{bamname}_flagstats.txt","w") as file:
        subprocess.run(cmd, stdout = file, stderr = True)

def pysamstats(m_df, min_base_q):
    for library in m_df.index: 
        pysamstats_helper(min_base_q, m_df, library)
    print("Pysamstats has completed for all libraries")

def pysamstats_helper(min_base_q, m_df, library):
    params = set_parameters(m_df, library)
    reference_dir = params.get("reference_dir")
    reference_fa = params.get("reference_fa")
    star_out_prefix = params.get("star_out_prefix") + "CHIMERIC_"
    working_dir = os.getcwd()
    cmd = ["docker" ,"run" , "--rm" ,
        "-v", f"{working_dir}:/data",
        "--workdir" , "/data",
        "pysamstats:1.3.0",
        "pysamstats" ,"--min-baseq",str(min_base_q), "--pad",
        "--max-depth=999999999", "-t", "variation_strand", 
        "-f", reference_dir + "/" + reference_fa, f"{star_out_prefix}Aligned.sortedByCoord.out.bam"]
    with open(f"{star_out_prefix}pysamstats.tsv","w") as file:
        subprocess.run(cmd, stdout = file, stderr = True)
    print("Pysamstats has completed for " + library)

def soft_clip(m_df, run, soft_clip_pct_t, omit_ends):
    """Function to get soft clipped positions to a df (soft_clip_count), calculate soft clipped pct at each position (soft_clip_table), plot soft clipped pct (soft_clip_plot)"""
    bam_list = []
    soft_clipped_width_list = []
    for library in m_df.index:
        bam_list.append(glob.glob(f"{library}/{library}*CHIMERIC_Aligned.sortedByCoord.out.bam")[0])
    for bamname in bam_list:
        if os.path.getsize(bamname) > 0:
            print (f"processing {bamname}")
            library = bamname.split('/')[0]
            bam_n = bamname.split('/')[1]
            bam_n = bam_n.split('_Aligned')[0]
            soft_clip_df = soft_clip_count(m_df, library, bamname)
            soft_clip_w_pct= soft_clip_table(library, run, bam_n, soft_clip_df)
            soft_clipped_width_df = soft_clip_plot(m_df, library, run, omit_ends, bam_n, soft_clip_w_pct, "pct")
            soft_clipped_width_list.append(soft_clipped_width_df)
    
    soft_clipped_dict = {}
    for l,df in zip(bam_list, soft_clipped_width_list):
        soft_clipped_dict[l] = df
    width_df = pd.concat(soft_clipped_dict)

    width_df_sorted = width_df.sort_values(by = ['library','peak pct'], ascending = [True, False])
    width_df_sorted = width_df_sorted[(width_df_sorted['peak pct'] > soft_clip_pct_t)&(width_df_sorted['coverage'] > 2500)]
    width_df_sorted.to_csv(f"results/soft_clipped_peak_df.csv") 
    return(width_df_sorted)


def soft_clip_count(m_df, library, bamname):
    """Function to get soft clipped positions to a df"""

    bam = pysam.AlignmentFile(bamname, 'rb')
    params = set_parameters(m_df,library)
    reference_dir = params.get("reference_dir")
    reference_fa = params.get("reference_fa")
    for record in SeqIO.parse(f"{reference_dir}/{reference_fa}", "fasta"):
        length = len(record.seq)

    # initiate an empty array with the same length as reference seq length
    arr=[]
    arr = [0 for i in range(length)] 

    for read in bam.fetch(until_eof=True):
        seq = read.query_sequence
        ref_pos = read.get_reference_positions(full_length = True)
        if len(ref_pos)>0:
            # if 5' is soft clipped
            if str(ref_pos[0]) == 'None':
                l_stack = []
                # use stack to track clipped position
                for pos in ref_pos:
                    if (str(pos) == 'None'):
                        l_stack.append(str(pos))
                    # when the next position is not clipped, end searching
                    elif (len(l_stack) != 0 and l_stack[-1] == 'None'):
                        l_stack.append(pos)
                        # breaks inner loop
                        break
                # store the info of soft clipped positions
                clipped_pos = l_stack.pop()
                if (isinstance(clipped_pos, int)):
                    # otherwise, the whole read is unmapped
                    while (len(l_stack) != 0):
                        clipped_pos -= 1
                        l_stack.pop()
                        if clipped_pos >= 0:
                            arr[clipped_pos] = arr[clipped_pos] + 1
                            
                
            #if 3' is soft clipped
            if str(ref_pos[-1]) == 'None':
                # counting from the reverse order
                ref_pos.reverse()
                r_stack = []
                for pos in ref_pos:
                    if (str(pos) == 'None'):
                        r_stack.append(str(pos))
                    elif (len(r_stack) != 0 and r_stack[-1] == 'None'):
                        r_stack.append(pos)
                        # breaks inner loop
                        break
                clipped_pos = r_stack.pop()
                if (isinstance(clipped_pos, int)):
                    while (len(r_stack) != 0):
                        clipped_pos += 1
                        r_stack.pop()
                        if clipped_pos < length:
                            arr[clipped_pos] = arr[clipped_pos] + 1
                            
    myRange = np.arange(0,len(arr),1)
    soft_clipped_df = pd.DataFrame({"position": myRange, 'soft_clipped_count':arr})
    return soft_clipped_df

def soft_clip_table(library, bam_n, soft_clip_df):
    columns = ["sample","ref","pos","coverage","nt","soft clipped count", "soft clipped pct"]
    soft_clipped_pct_df = pd.DataFrame(columns = columns)
    # read in pysamstats data
    file = glob.glob("{0}/{0}*CHIMERIC_pysamstats.tsv".format(library))
    pysam_df = pd.read_csv(file[0],sep ="\t")
    i = 0
    while i < len(pysam_df):
        # get necessary values
        # pysamstats pos starts from 1 while soft_clip_df pos starts from 0
        position = i+1
        reads = max(pysam_df["reads_all"].iloc[i],1)
        ref = pysam_df["ref"].iloc[i]
        soft_clipped_count = soft_clip_df["soft_clipped_count"].iloc[i]
        soft_clipped_pct = round(soft_clipped_count / (reads + soft_clipped_count) * 100,2)
        sample = library 
        sq = pysam_df["chrom"].iloc[i]
        soft_clipped_pct_df.loc[len(soft_clipped_pct_df.index)] = [sample,sq,position,reads,ref,
                                            soft_clipped_count, soft_clipped_pct]
        i+=1

    soft_clipped_pct_df.to_csv(f"results/{bam_n}_soft_clipped_pct.csv", index = False)
    print (f"results/{bam_n}_soft_clipped_pct.csv")
    return soft_clipped_pct_df
        
def soft_clip_plot(m_df, library, run, omit_ends, bam_n, soft_clip_w_pct, display):
    #read in pysamstats data
    file = glob.glob("{0}/{0}*CHIMERIC_pysamstats.tsv".format(library))
    pysam_df = pd.read_csv(file[0],sep ="\t")
    ref = m_df.loc[m_df.index==library,"Expected Reference"][0]

    #plot
    sns.set_style('white')
    sns.set_context('poster')
    fig = plt.figure(figsize=[16,10])

    ax1 = fig.add_subplot(111)
    start = 50
    end = 150
    #As of now leave as reads_all
    pysam_df.iloc[int(start):len(pysam_df) - int(end)]['reads_all'].plot()
    #Coverage plot
    ax1.get_yaxis().get_major_formatter().set_scientific(False)
    ax1.set_ylim(0)
    ax1.set_ylabel("Coverage ",{"size":24})
    ax1.set_xlabel(f"Nucleotide Number {int(start)} to {len(pysam_df) -int(end)}", {"size":24}) if omit_ends else ax1.set_xlabel("Nucleotide Number", {"size":24})
    ax1.set_title(f"{bam_n} {ref}"  ,{"size":30})

    #Soft clipped plot
    ax2= ax1.twinx()
    ax2.set_ylabel("% Soft Clipped", {"size": 24, "rotation": 270},labelpad =20)
    ax2.set_ylim(0,30)
    ax2.set_xlim(start,len(pysam_df)-end)

    soft_clipped_df = soft_clip_w_pct['soft clipped pct']
    soft_clipped_df.iloc[int(start):len(pysam_df)-int(end)].plot(color= 'red')

    #smooth curve
    x_avg = mean(soft_clipped_df[int(start):(len(pysam_df) -int(end))])
    xhat = savgol_filter(soft_clipped_df, 51, 3)
    xhat_ser = pd.Series(xhat, copy=False) 
    #calling peak and width
    peaks, _ = find_peaks(xhat_ser[int(start):(len(pysam_df) -int(end))], prominence=x_avg/2)
    results_half = peak_widths(xhat_ser[int(start):(len(pysam_df) -int(end))], peaks, rel_height=0.5)
    # add back int(start)
    for i in range (0,len(peaks)):
        peaks[i] = int(peaks[i]) + int(start)
    for j in range(1,len(results_half[1:])):
        for m in range (0, len(results_half[1:][j])):
            results_half[1:][j][m] = int(results_half[1:][j][m]) + int(start)
    for n in range(0,len(results_half[0])): 
        results_half[0][n] = int (results_half[0][n])
    # store width into df
    values = []
    coverage = []
    for p in peaks:
        # peak is not necessarily the highest point, so look at nearby pos
        values.append(max(soft_clip_w_pct.loc[p-2,'soft clipped pct'], soft_clip_w_pct.loc[p-1,'soft clipped pct'], soft_clip_w_pct.loc[p,'soft clipped pct'], soft_clip_w_pct.loc[p+1,'soft clipped pct'], soft_clip_w_pct.loc[p+2,'soft clipped pct']))
        coverage.append(soft_clip_w_pct.loc[p,'coverage'])
    soft_clipped_width_df = pd.DataFrame()
    soft_clipped_width_df["library"] = [library]*len(peaks)
    soft_clipped_width_df["start"] = results_half[1:][1]
    soft_clipped_width_df["end"] = results_half[1:][2]
    soft_clipped_width_df["width"] = results_half[0]
    soft_clipped_width_df["soft clipped peak pos"] = peaks
    soft_clipped_width_df["coverage"] = coverage
    soft_clipped_width_df["peak pct"] = values
    sort_soft_clipped_width_df = soft_clipped_width_df.sort_values(by = 'peak pct', ascending=False)
    

    # plot elements
    xhat_ser[int(start):(len(pysam_df) -int(end))].plot(color='orange')
    plt.plot(peaks, xhat_ser[peaks], ".",markersize=30,color='orange')
    plt.hlines(*results_half[1:], color="C2")
    ax2.get_yaxis().get_major_formatter().set_scientific(False)

    #legend
    cov_patch = mpatches.Patch(color = 'blue', label='Coverage')  
    soft_clipped_patch = mpatches.Patch(color='red', label='Soft Clipped Percentage')
    soft_clipped_marker_patch = mpatches.Patch(color='orange', label='Soft Clipped Smoothed and Peak')
    soft_clipped_width_patch = mpatches.Patch(color="C2", label='Soft Clipped Width')

    lgd = fig.legend(
    handles=[cov_patch,soft_clipped_patch,soft_clipped_marker_patch,soft_clipped_width_patch],
    title="",
    bbox_to_anchor=(1.3,1),
    loc="upper right")
    
    fig.savefig(f"results/{library}_soft_clipped_{display}.png",bbox_inches='tight')
    return (sort_soft_clipped_width_df)


def generate_report(run, s3_suffix, omit_ends, report_threshold):
    end_time = time.monotonic()
    s3_dir = f"{run}_{s3_suffix}" if s3_suffix else f"{run}"
    job_name = f"{run}"
    omit_param = "TRUE" if omit_ends else "FALSE"
    cmd = ["Rscript", "-e", f"rmarkdown::render('SoftClip_report.Rmd', output_file = 'SoftClip_report.html', params= list(job=\"{job_name}\",s3_dir=\"{s3_dir}\",omit ={omit_param}, report_threshold = {report_threshold}))"]
    with open("report.log", "w") as f:
        print(subprocess.run(cmd,stdout = True, stderr = f))

def object_to_s3_url(object_name, bucket_string):
    s3_url = f'https://{bucket_string}/{object_name}'
    return s3_url  

def setup_local_backup(send_list, rename_list, send_s3_dir):
    os.makedirs(send_s3_dir, exist_ok = True)
    for s,r in zip (send_list, rename_list):
        shutil.copy(s, f'{send_s3_dir}/{r}')

def get_s3_urls(rename_list, send_s3_dir, bucket_string, s3_output_prefix):
    url_list = [f'https://{bucket_string}/{s3_output_prefix}/{send_s3_dir}/{r}' for r in rename_list]
    return(url_list)
    

def sync_to_s3(send_s3_dir, bucket_string, s3_output_prefix):
    aws_command = f'aws s3 sync {send_s3_dir} s3://{bucket_string}/{s3_output_prefix}/{send_s3_dir}/'
    proc = subprocess.run(aws_command, shell=True)

def dir_to_s3(send_list, rename_list, send_s3_dir, bucket_string, s3_output_prefix):
    setup_local_backup(send_list, rename_list, send_s3_dir)
    sync_to_s3(send_s3_dir, bucket_string, s3_output_prefix)
    url_list = get_s3_urls(rename_list, send_s3_dir, bucket_string, s3_output_prefix)

def backup_s3(bucket_string, s3_output_prefix):
    # soft_clip
    send_s3_dir = "soft_clip"
    send_list = glob.glob(f"results/*soft_clipped_pct.csv")
    rename_list = [os.path.basename(f) for f in send_list]
    dir_to_s3(send_list, rename_list, send_s3_dir, bucket_string, s3_output_prefix)

    #reference files
    sync_to_s3("fasta", bucket_string, s3_output_prefix)

    #summary
    send_s3_dir = "summary"
    os.makedirs(send_s3_dir, exist_ok = True)
    shutil.copy("results/soft_clipped_peak_df.csv", f'{send_s3_dir}/')
    shutil.copy("SoftClip_report.html", f'{send_s3_dir}/')
    sync_to_s3(send_s3_dir, bucket_string, s3_output_prefix)