import pandas
import os
import sys

#different readlengths for trimmomatic outputs
READLENGTHS = [36,75]

#from utils import metautils
#metautils.srpMeta('SRP142360')
#study_accession            experiment_accession            sample_accession            run_accession            experiment_title            experiment_attribute            taxon_id            library_selection            library_layout            library_strategy            library_source            library_name            bases            spots            adapter_spec            avg_read_length
#      SRP141397                      SRX3980107                  SRS3205258                SRR7049033                S78_shotgun                                                   0                       RANDOM                  PAIRED -                         WGS               METAGENOMIC             S78_shotgun        229716627          1009555                                 227.54245880610765

	
#study_accession    experiment_accession                experiment_title  experiment_desc   organism_taxid    organism_name     library_strategy  library_source    library_selection                   sample_accession  sample_title      instrument        total_spots       total_size        run_accession     run_total_spots   run_total_bases   run_alias         sra_url_alt1      sra_url_alt2      sra_url           experiment_alias  isolation_source  collection_date   geo_loc_name      lat_lon           ref_biomaterial   BioSampleModel    host              ena_fastq_url_1   ena_fastq_url_2   ena_fastq_ftp_1   ena_fastq_ftp_2   full1             full2             pair1             pair2             
#SRP141397          SRX3980708                          MP74_16S          MP74_16S          256318            metagenome        AMPLICON          METAGENOMIC       PCR                                 SRS3205645        N/A               Illumina MiSeq    550               257302            SRR7049034        550               276100            MP74_2.fastq.gz   https://storage.googleapis.com/sra-pub-src-4/SRR7049034/MP74_2.fastq.gz https://sra-pub-src-4.s3.amazonaws.com/SRR7049034/MP74_2.fastq.gz       https://sra-downloadb.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-14/SRR7049034/SRR7049034.1  N/A               N/A               25-Oct-2017       USA:Philadelphia  39.95 N 75.16 W   MP74              Metagenome or environmental         Homo sapiens      http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR704/004/SRR7049034/SRR7049034_1.fastq.gz           http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR704/004/SRR7049034/SRR7049034_2.fastq.gz           era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR704/004/SRR7049034/SRR7049034_1.fastq.gz        era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR704/004/SRR7049034/SRR7049034_2.fastq.gz        raw/SRP141397/SRX3980708/SRR7049034_1.fastq           raw/SRP141397/SRX3980708/SRR7049034_2.fastq           sratofastq/SRR7049034_1.fastq.gz    sratofastq/SRR7049034_2.fastq.gz

#read in an SPR metadata

class srpMeta():
    def __init__(self,SRP):
        self.srp=SRP
        self.st = pandas.read_csv("metadata/"+SRP+".metadata",sep="\t")
        #cell.line not pandas friendly
        self.st.columns = self.st.columns.str.replace('.', '')
        #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz', 
    
    def process(self,gstype='gs'):
        #create the gs specific paths
        #gs://truwl-quertermous/SRR7058289/102901.1.fastq.gz
        #gs://truwl-quertermous/SRR7058289/102901.2.fastq.gz
        from utils.metautils import srpMeta
        srp=srpMeta('SRP142360').process()
        if gstype=='https':
            prefix="https://storage.googleapis.com/truwl-quertermous"
        elif gstype=='gs':
            prefix="gs://truwl-quertermous/"
        else:
            raise ValueError("need a gstype of https or gs")
        self.st["truwl_pair1"] = "{0}/{1}/{2}.1.fastq.gz".format(prefix,self.st['run_accession'],self.st['cell.line'])
        self.st["truwl_pair2"] = "{0}/{1}/{2}.2.fastq.gz".format(prefix,self.st['run_accession'],self.st['cell.line'])

    def getSRRs(self):
        return(self.st['run_accession'].tolist())
    
    def getATACFiles(self,flatten=True,compress=True):
        """
        flatten - all files go in top directory vs SRP/EXP/SRR
        compress - .gz suffix
        """
        files=[]
        filesofinterest=self.st[self.st['library_strategy']=='ATAC-seq'].copy()
        files=list(filesofinterest['truwl_pair1'].astype(str))+list(+filesofinterest['truwl_pair2'].astype(str))
        
    def getRNASEQFiles(self,flatten=True,compress=True):
        """
        flatten - all files go in top directory vs SRP/EXP/SRR
        compress - .gz suffix
        """
        files=[]
        filesofinterest=self.st[self.st['library_strategy']=='ATAC-seq'].copy()

    def getATACManifest(self):
        cell_lines = self.st['cell.line'].tolist()
#         {
#     "atac.pipeline_type" : "atac",
#     "atac.genome_tsv" : "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv",
#     "atac.fastqs_rep1_R1" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF341MYG.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF106QGY.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep1_R2" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF248EJF.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF368TYI.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep2_R1" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF641SFZ.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF751XTV.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF927LSG.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF859BDM.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF193RRC.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF366DFI.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep2_R2" : [
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF031ARQ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF590SYZ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF734PEQ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF007USV.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF886FSC.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF573UXK.subsampled.400.fastq.gz"
#     ],
#     "atac.paired_end" : true,
#     "atac.auto_detect_adapter" : true,
#     "atac.enable_xcor" : true,
#     "atac.title" : "ENCSR356KRQ (subsampled 1/400)",
#     "atac.description" : "ATAC-seq on primary keratinocytes in day 0.0 of differentiation"
# }
    
    def getAllFiles(self):
        getFilesFromRunList(self.populationRuns('ALL'))

    #>>> metautils.calculateExpected('ALL')
    #Expecting 413 fastq and 253 bams
    def calculateExpected(self,population):
        fastq = self.getFilesFromRunList(self.populationRuns(population))
        bams = self.populationRuns(population)
        print("Expecting {0} fastq and {1} bams".format(len(fastq),len(bams)))
    
    #GBR FIN ALL YRI TSI 
    #>>> metautils.populationRuns('ALL')
    #['NA12286.2.MI_120126_5', 'NA10851.4.M_120208_1',
    def populationRuns(self,population):
        return(self.getSRRs())

    #Check to see if paired-end sequencing was used in each sample so appropriate stAR parameters are called
    def isPaired(self,run):
        if self.st.loc[self.st['Assay Name']==run]['Comment[LIBRARY_LAYOUT]'].all()=='PAIRED':
            return(True)
    
    #Either run star with paired or unpaired options
    def starProgram(self,run):
        if self.isPaired(run):
            return("star_paired.sh")
        else:
            return("star_unpaired.sh")
    
    #Run trimmomatic bash script with paired or unpaired options
    def trimProgram(self,run):
        if self.isPaired(run):
            return("trimmomatic_paired.sh")
        else:
            return("trimmomatic_unpaired.sh")
            
    #pairedOrSinglefastqInput('NA06984.1.M_111124_4')
    #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz']
    #pairedOrSinglefastqInput('NA06985.1.MI_120104_3')
    #['NA06985.1.MI_120104_3_1.fastq.gz']
    
    #set file names for trimmomatic based on whether a trimmed input or output is expected
    def pairedOrSinglefastqInput(self,run,trimmed=False,readlength=None):
        if trimmed:
            extName=".trimmed_{0}.fastq.gz".format(readlength)
        else:
            extName=".fastq.gz"
        if self.isPaired(run):
            return([run+"_1"+extName,run+"_2"+extName])
        else:
            return([run+"_1"+extName])
    

    def getPathFromName(self,name):
        return(ALL_LOOKUP[name])
    
    #ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188404/ERR188404_1.fastq.gz

    
    def getType(self,runName):
        return(self.getProperty(runName, 'Platform'))
    
    
    def getProperty(self, runName, column):
        stseries = self.st.loc[self.st['Run'] == runName, column]
        if len(self.stseries) == 1:
            return(self.stseries.to_self.string(index=False))  # pandas is ridiculous
        elif len(self.stseries) > 1:
            raise ValueError("Not expecting to match multiple run names, just 1")
        else:
            raise ValueError("Can't find that run {0}".format(runName))
    
    
    def getMemory(self, runName):
        return(int(getProperty(runName, 'size_MB')))
    
    #Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed
    def cleanJobname(self, name):
        return(name.replace('.','_'))
    
    def getECS(self, runName, units, program):
        if program == 'minimap':
            mb = 32000
        elif program == 'star':
            superIntensive = ['NA12716.7.M_120219_6', 'NA12775.7.M_120219_5',
                              'NA12814.7.M_120219_3', 'NA11831.7.M_120219_2',
                              'NA11994.2.M_111216_7', 'NA11893.7.M_120219_3']
            if runName in superIntensive:
                mb = 64000
            else:
                mb = 48000
        elif program == 'IsoModule':
            mb = 192000
        elif program == 'samtoolsindex':
            mb = 16000
        elif program == 'samtoolsmerge':
            mb = 16000
        elif program == 'samtoolssubsample':
            mb = 16000
        elif program == 'fastx':
            mb = 8000
        elif program == 'trimmomatic':
            mb = 8000
        elif program == 'rmats':
            mb = 16000
        else:
            raise ValueError
        if units == 'bytes':
            return(1048576*mb)
        elif units == 'mb':
            return(mb)
    # 2  # number of sample
    # 3  # number of replicates in sample #1
    # sam1_rep1.bam
    # sam1_rep2.bam
    # sam1_rep3.bam
    # 3  # number of replicates in sample #2
    # sam2_rep1.bam
    # sam2_rep2.bam
    # sam2_rep3.bam
    
    #Used for a two way comparison between various samples
    def twoSampleComparisonManifest(self, biosamp1, biosamp2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for biosamp in [biosamp1, biosamp2]:
            run = getRunBamsFromBioSample(biosamp, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    
        #Used for a two way comparison between various samples
    def twoTreatmentComparisonManifest(self, treatment1, treatment2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for treatment in [treatment1, treatment2]:
            run = getRunBamsFromTreatment(treatment, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    # "panorama-clk-repro/SRP091981/
    
    #>>> metautils.getRunBamsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam.bai']
    #Get a List of bams originating from the same biological source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunBamsFromBioSample(self, biosamp, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromBioSample(biosamp)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    
    def getRunBamsFromTreatment(self, treatment, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromTreatment(treatment)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    #>>> metautils.getRunsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4', 'NA12778.1.M_111124_4', 'NA12778.1.MI_120104_3']
    #Get a List of runs from the same source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunsFromBioSample(self, biosample,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['Source Name'] == biosample) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))
    
    def getRunsFromTreatment(self, treatment,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['treatment'] == treatment) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))

    #metautils.populationBiosamples('ALL')
    #['NA12778', 'NA12045', 'NA12144', ...
    #Returns either all runs that pass the QC check for a specific population or all populations depending on arguments
    def populationBiosamples(self,population):
        if(population=='ALL'):
            return(list(set(self.st.loc[(self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
        else:
            return(list(set(self.st.loc[(self.st['Characteristics[population]'] == population) & (self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
    
    #This returns all sequencing runs over 70 nt (longruns) looks through the population being considered appends them to the sample List using the mergedbam file naming convention, and then appends all runs from that population with a 36nt trimmed suffix
    def getListOfLengthOutputs(self,population):
        x = 0
        biosamps = self.populationBiosamples(population)
        biosamp_lengthList = []
        longruns = list(set(self.st.loc[self.st['Comment[SEQUENCE_LENGTH]'].astype(int)<70]['Source Name'].toList()))
        for item in longruns:
            if (item in biosamps) and (item != 'NA07000'):
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        for item in biosamps:
            if item != 'NA07000':
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        return(biosamp_lengthList)

if __name__== "__main__":
    srpMeta(sys.argv[1]).get16SFiles()
