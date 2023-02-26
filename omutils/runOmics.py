import argparse
import subprocess
import glob
import os
from collections import namedtuple


class Omics:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description="Processing multi-sample QC and quanlification for raw read of RNA-seq.")
        parser.add_argument("-i", "--input", type=str,
                            help="Specify your input directory including raw fastq files")
        parser.add_argument("-f", "--suffix1", type=str, default=".r1",
                            help="Specify suffix of read1, such as '.r1' for sample1_r1.fq. <PE model only>")
        parser.add_argument("-b", "--suffix2", type=str, default=".r2",
                            help="Specfiy suffix of read2, such as '.r2' for sample1_r2.fq. <PE model only>")
        parser.add_argument("-o", "--output", type=str, default=".",
                            help="Sepecfiy output directory. <default = ./>")
        parser.add_argument("-g", "--genome", type=str,
                            help="Path to genome fasta")
        parser.add_argument("-gt", "--gtf", type=str,
                            help="Path to genome gtf file")
        parser.add_argument("-gi", "--indexdir", type=str,
                            help="Specify a directory to create/include genome index")
        parser.add_argument("-t", "--thread", type=int, default=1,
                            help="Specify thread you want to use. <default = 1>")
        parser.add_argument("-z", action="store_true",
                            help="Specify your fq files are compressed")
        parser.add_argument("-p", default=False, action="store_true",
                            help="Specify read is paired")
        parser.add_argument("-s", default=False, action="store_true",
                            help="Specfiy read is single")
        parser.add_argument(
            "-om", "--omics", default="rnaseq", help="Spefify omics")
        parser.add_argument("-m", "--module", default="qc", help="Run module")
        parser.add_argument("-idb", "--genomeSAindexNbases", type=int, default=14, help="Specify genomeSAindexNbases for STAR to create index")
        self.args = parser.parse_args()
        self.input_d = self.args.input
        self.out_d = self.args.output
        self.thread = self.args.thread
        self.genome = self.args.genome
        self.gtf = self.args.gtf
        self.indexdir = self.args.indexdir
        self.suffix1 = self.args.suffix1
        self.suffix2 = self.args.suffix2
        self.omics = self.args.omics
        self.module = self.args.module
        self.qcdir = os.path.join(self.out_d, "QC")
        self._out = ["are", "is"]
        self.outcounts = os.path.join(self.out_d, "featureCounts.txt")
        self.aligndir = os.path.join(self.out_d, "aligndir")
        self.countsdir = os.path.join(self.out_d, "countsdir")
        self.genomeSAindexNbases = self.args.genomeSAindexNbases
        self.qcfilelist = []
        print(
            f"Your input directory is {self.input_d}, output directory is {self.out_d}, \
            thread {self._out[0] if self.thread >1 else self._out[1]} {self.thread}.")
        if not self.input_d:
            raise TypeError("Please specify input directory!")
        if self.args.s and (self.suffix1 != ".r1" or self.suffix2 != ".r2"):
            raise TypeError("Only PE model needs to specify -f and -b")
        if self.args.p and not self.args.s:
            print("Model: paired read")
        elif self.args.s and not self.args.p:
            print("Model: single read")
        elif self.args.p and self.args.s:
            raise TypeError("Please specify your read type!")
        if not os.path.isdir(self.out_d):
            os.mkdir(self.out_d)

    def runQC(self):
        if not os.path.isdir(self.qcdir):
            os.mkdir(self.qcdir)
        subprocess.Popen(
            "for id in %s/*fastq*; do mv $id ${id/fastq/fq}; done" % self.input_d, shell=True)
        for sample in glob.glob(f"{self.input_d}/*.fq*"):
            gzsuffix = ".fq.gz" if sample.endswith(".fq.gz") else ".fq"
            if self.args.s:
                self._qcnamedt = namedtuple("qcnamedt", ["r1in", "r1out"])
                sample_name = os.path.basename(sample).split(f".fq.gz")[0] if sample.endswith(
                    ".fq.gz") else os.path.basename(sample).split(".fq")[0]
                outfile = os.path.join(
                    self.qcdir, sample_name+".clean"+gzsuffix)
                self.qcfilelist.append(self._qcnamedt(sample, outfile))
                print(f"Sample file is {sample}, outfile is {outfile}")
                # subprocess.Popen(f"fastp -i {sample} -o {outfile}", shell=True)
                os.system(f"fastp -i {sample} -o {outfile} -j {self.qcdir}/{sample_name}.json -h {self.qcdir}/{sample_name}.html")
            else:
                if self.suffix1 in sample:
                    self._qcnamedt = namedtuple(
                        "qcnamedt", ["r1in", "r1out", "r2in", "r2out"])
                    sample_name = os.path.basename(sample).split(f"{self.suffix1}.fq.gz")[0] if sample.endswith(
                        ".fq.gz") else os.path.basename(sample).split(f"{self.suffix1}.fq")[0]
                    read1_outfile = os.path.join(
                        self.qcdir, sample_name+self.suffix1+".clean"+gzsuffix)
                    read2_infile = os.path.join(os.path.dirname(
                        sample), sample_name+self.suffix2+gzsuffix)
                    read2_outfile = os.path.join(
                        self.qcdir, sample_name+self.suffix2+".clean"+gzsuffix)
                    self.qcfilelist.append(self._qcnamedt(
                        sample, read1_outfile, read2_infile, read2_outfile))
                    print(
                        f"read1 file is {sample}, outfile is {read1_outfile}, read2 file is {read2_infile}, \
                            outfile is {read2_outfile}")
                    # subprocess.Popen(f"fastp -i {sample} -o {read1_outfile} -I {read2_infile} -O {read2_outfile}", shell=True)
                    os.system(
                        f"fastp -i {sample} -o {read1_outfile} -I {read2_infile} -O {read2_outfile} -j {self.qcdir}/{sample_name}.json -h {self.qcdir}/{sample_name}.html")

    def buildIndex(self):
        if not os.path.isdir(self.indexdir):
            os.mkdir(self.indexdir)
        os.system(f"STAR --genomeSAindexNbases {self.genomeSAindexNbases} --runMode genomeGenerate --genomeDir {self.indexdir} \
                  --runThreadN {self.thread} --genomeFastaFiles {self.genome} --sjdbGTFfile {self.gtf}")

    def align(self):
        if not os.path.isdir(self.aligndir):
            os.mkdir(self.aligndir)
        for qcfile in self.qcfilelist:
            samplename = os.path.basename(qcfile.r1out.split(f"{self.suffix1}")[0])
            alignoutfile = os.path.join(self.aligndir, samplename)
            print(qcfile)
            if self.args.p:
                rfin = qcfile.r1out + " " + qcfile.r2out
            else:
                rfin = qcfile.r1out
            os.system(f"STAR --genomeDir {self.indexdir} --runThreadN {self.thread} \
                        --readFilesIn {rfin} \
                        --readFilesCommand zcat \
                        --outFileNamePrefix {alignoutfile} \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMstrandField intronMotif \
                        --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical \
                        --outBAMsortingThreadN {self.thread}")

    def quantify(self):
        if not os.path.isdir(self.countsdir):
            os.mkdir(self.countsdir)
        bamfiles = " ".join(glob.glob(f"{self.aligndir}/*.bam"))
        os.system(f"featureCounts -t exon -g gene_id \
                    -Q 10 --primary -s 0 -p -T {self.thread} \
                    -a {self.gtf} \
                    -o {self.countsdir}/featureCounts.txt \
                    {bamfiles}")

    def pipeline(self):
        self.runQC()
        self.buildIndex()
        self.align()
        self.quantify()


if __name__ == "__main__":
    omics = Omics()
    omics.pipeline()
