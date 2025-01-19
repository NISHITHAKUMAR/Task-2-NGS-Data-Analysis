##Task1 Python script
import os 
import subprocess
import shutil
import os
from urllib.parse import urlparse
from collections import defaultdict

def CREATE_FILE_PAIR(DATA_DIR):
    result_dict = defaultdict(lambda: ["", ""])
    for path in os.listdir(DATA_DIR):
        base_name = path.replace('_R1_001', '').replace('_R2_001', '')
        if '_R1_' in path:
            result_dict[base_name][0] = path
        elif '_R2_' in path:
            result_dict[base_name][1] = path
    result_dict = dict(result_dict)
    return result_dict

def run_subprocess(command):
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        rc = process.wait()
        if rc != 0:
            print("Process exited with non-zero status")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

def execute_fastaq(pair1_entry, pair2_entry, thr):
        qual_path = os.path.join(os.getcwd(), "Quality_Control")
        os.makedirs(qual_path, exist_ok=True)
        

        if os.path.isfile(pair1_entry) and os.path.isfile(pair2_entry):
            command = ["fastqc", pair1_entry, pair2_entry, "-t", str(int(thr)), "-o", qual_path]
            run_subprocess(command)
        else:
            print("Path Error", "Please select Path properly")

        CMD_FASTP = f"""fastp -i {pair1_entry} -w 10 -o {os.path.join(qual_path, f"trim_{os.path.basename(pair1_entry)}")} -I {pair2_entry} -O {os.path.join(qual_path, f"trim_{os.path.basename(pair2_entry)}")} -q 30"""
        os.system(CMD_FASTP)

        command_2 = ["fastqc", os.path.join(qual_path, f"trim_{os.path.basename(pair1_entry)}"), os.path.join(qual_path, f"trim_{os.path.basename(pair2_entry)}"), "-t", str(int(thr)), "-o", qual_path]
        
        print(command_2)
        run_subprocess(command_2)



def download_genome(genome_entry):
    os.makedirs(os.path.join(os.getcwd(), "Ref_Genome"), exist_ok=True)
    os.chdir(os.path.join(os.getcwd(), "Ref_Genome"))
    if os.path.isfile(genome_entry):
        return(genome_entry)
    else:
        name = genome_entry.split("/")[-1].strip()
        parsed = urlparse(genome_entry)
        if all([parsed.scheme, parsed.netloc]):
            os.system(f"curl {genome_entry} -o {name}")
            os.system(f"gunzip {name}")
        return os.path.join(os.getcwd(), name.replace(".gz", ""))


def build_Genome_Indexing(fa_path):

    try:
        os.chdir(os.path.dirname(fa_path))
    except:
        pass

    name = fa_path.split("/")[-1].strip()
    index_files = os.listdir(os.path.dirname(fa_path))  
    if len(index_files) >= 2:
        print("Indexing was already done.")
        return fa_path
    else:
        print("Buildin genome Index ..............")
        # CMD = ['bwa', 'index',  fa_path]
        os.system(f"bwa index {fa_path}")
        return fa_path

def BWA_MEME(pair1_entry, pair2_extry, fa_path):

    common_filename = os.path.commonprefix([os.path.basename(pair1_entry), os.path.basename(pair2_extry)])
    CMD = f"bwa mem {fa_path} {pair1_entry} {pair2_extry} >{common_filename}.sam"
    os.system(CMD)
    return(os.path.join(os.getcwd(), f"{common_filename}.sam"))

def execute_MUTATE2(PATH, TUMOR_BAM, NORMAL_BAM):
    REF_VCF_DIR = os.path.join(PATH, 'Ref_VCF')
    MUTACT_2_PATH = os.path.join(PATH, "NGS_Tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar Mutect2")
    NORMAL_NAME  = os.path.splitext(os.path.basename(NORMAL_BAM))[0]
    REF_GENOME = os.path.join(PATH, "Ref_Genome/genome.fa")

    CMD  = f"/usr/lib/jvm/java-1.21.0-openjdk-amd64/bin/java -jar {MUTACT_2_PATH} -R {REF_GENOME} -I {TUMOR_BAM} -I {NORMAL_BAM} -normal {NORMAL_NAME} --germline-resource {os.path.join(REF_VCF_DIR, 'af-only-gnomad.vcf.gz')} --panel-of-normals {os.path.join(REF_VCF_DIR, 'GermlineHetPon.38.vcf.gz')} -O mutect2_somatic.vcf"
    try:
        os.system(CMD)
    except Exception as e:
        print(e)

    CMD_2 = f"/usr/lib/jvm/java-1.21.0-openjdk-amd64/bin/java -jar {MUTACT_2_PATH} -R {REF_GENOME} -I {NORMAL_BAM} -O mutect2_somatic.Normal.vcf"
    try:
        os.system(CMD_2)
    except Exception as e:
        print(e)


added_file_paths = []
def process_sam_file(file_path):
    results_dir = os.path.dirname(file_path)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    bam_file = os.path.join(results_dir, f"{base_name}.bam")
    markdup_bam_file = os.path.join(results_dir, f"{base_name}.markdup.bam")
    sample_added_bam = os.path.join(results_dir, f"{base_name}.AddSample.markdup.bam")

    try:
        # Convert SAM to BAM and sort
        subprocess.run(f"samtools view -bS {file_path} | samtools sort --threads 10 -o {bam_file}", shell=True, check=True)
        print(f"Converted and sorted {file_path} to {bam_file}")

        # Mark duplicates
        # subprocess.run(f"samtools markdup {bam_file} {markdup_bam_file}", shell=True, check=True)
        os.system(f"picard-tools MarkDuplicates I={bam_file} O={markdup_bam_file} M=marked_dup_metrics.txt")
        print(f"Marked duplicates in {bam_file} to {markdup_bam_file}")

        CMD = f"picard-tools AddOrReplaceReadGroups I={markdup_bam_file} O={sample_added_bam} RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 SM={base_name}"
        os.system(CMD)
        # Index the marked duplicates BAM file
        subprocess.run(f"samtools index {sample_added_bam}", shell=True, check=True)
        print(f"Indexed {sample_added_bam}")

        added_file_paths.append(sample_added_bam)

        # ## vatriation calling using mutete2 
        # CMD_2 = f"/usr/lib/jvm/java-1.21.0-openjdk-amd64/bin/java -jar {os.path.join(PATH, 'NGS_Tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar Mutect2')} -R {os.path.join(PATH, 'Ref_Genome/genome.fa')} -I PA220KH-lib09-P19-Tumor_S2_L001_R.markdup.bam --germline-resource {os.path.join(PATH,'Ref_VCF/af-only-gnomad.hg38.vcf.gz')} --panel-of-normals {os.path.join(PATH,'Ref_VCF/GermlineHetPon.38.vcf.gz')} -O {base_name}_.vcf"
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing the file: {e}")


# Set the working path
PATH = os.getcwd()
os.chdir(PATH)

# Directory containing raw data
DATA_DIR = os.path.join(PATH, 'Raw_DATA')

# Assuming CREATE_FILE_PAIR is a function that generates a dictionary of file pairs
FILE_DICT = CREATE_FILE_PAIR(DATA_DIR)
# Genome reference URL
#genome = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz"

# Download the genome reference
# try:
#     fapath = download_genome(genome)
# except Exception as e:
#     print(e)

# # Build genome indexing
# try:
#     build_Genome_Indexing(fapath)
# except Exception as e:
#     print(e)




fapath = "Ref_Genome/genome.fa"
os.makedirs(os.path.join(PATH, "Results"),exist_ok=True )
os.chdir(os.path.join(PATH, "Results"))

result_path = os.path.join(PATH, "Results")


# Loop through each key in FILE_DICT
for key in FILE_DICT:
    try:
        # Get the pair of files for the current key
        pair1_entry = os.path.join(DATA_DIR, FILE_DICT[key][0])
        pair2_entry = os.path.join(DATA_DIR, FILE_DICT[key][1])

        # Print the file paths for verification
        print(f"Processing: {key}")
        print(f"Pair 1: {pair1_entry}")
        print(f"Pair 2: {pair2_entry}")
        
        # Execute the fastq processing function
        execute_fastaq(pair1_entry, pair2_entry, 10)

        pair_1_trimmed = os.path.join(result_path, f"Quality_Control/trim_{os.path.basename(pair1_entry)}")
        pair_2_trimmed = os.path.join(result_path, f"Quality_Control/trim_{os.path.basename(pair2_entry)}")

        # Align with BWA MEME
        SAM_PATH = BWA_MEME(pair_1_trimmed, pair_2_trimmed, fapath)

        # Process the resulting SAM file
        process_sam_file(SAM_PATH)

    except Exception as e:
        print(f"Error processing {key}: {e}")

TUMOR_BAM_ADDED = [f for f in added_file_paths if 'Tumor' in os.path.basename(f)][0]
NORMAL_BAM_ADDED = [f for f in added_file_paths if 'Norm' in os.path.basename(f)][0]

try:
    execute_MUTATE2(PATH, TUMOR_BAM_ADDED, NORMAL_BAM_ADDED)
except Exception as e:
    print(e)



