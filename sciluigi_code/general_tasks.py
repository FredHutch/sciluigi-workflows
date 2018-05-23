import os
import sciluigi as sl

class LoadFile(sl.ExternalTask):
    path = sl.Parameter()

    def out_file(self):
        return sl.ContainerTargetInfo(self, self.path)


class FastqpTask(sl.ContainerTask):

    # Input: FASTQ
    in_fastq = None

    # Output is a summary of the FASTQ quality
    summary_path = sl.Parameter()

    input_mount_point = sl.Parameter(default="/mnt/input/")
    output_mount_point = sl.Parameter(default="/mnt/output/")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/fastqp:fastqp-v0.2"

    def out_summary(self):
        return sl.ContainerTargetInfo(self, self.summary_path)

    def run(self):

        input_targets = {
            "fastq": self.in_fastq()
        }

        output_targets = {
            "summary_file": self.out_summary()
        }

        self.ex(
            command="fastqp " +
                    "-e $summary_file " + 
                    "$fastq",
            input_targets=input_targets,
            output_targets=output_targets,
            input_mount_point=self.input_mount_point,
            output_mount_point=self.output_mount_point,
        )


class AlignFastqTask(sl.ContainerTask):

    # Inputs: FASTQ and reference database
    in_fastq = None
    in_ref_fasta = None

    # Parameter: Short name for reference
    ref_name = sl.Parameter()
    
    # Parameter: AWS S3 folder for this project
    base_s3_folder = sl.Parameter()

    # Parameter: Name for output file(s)
    sample_name = sl.Parameter()

    # Parameter: Number of threads for alignment
    threads = sl.Parameter()

    # Parameter: Temporary folder to use on the device
    temp_folder = sl.Parameter()

    # URL of the container
    container = "quay.io/fhcrc-microbiome/bwa:v0.7.17--3"

    def out_bam(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.base_s3_folder,
                "align_bwa_{}".format(self.ref_name),
                "{}.{}.bam".format(self.ref_name, self.sample_name)
            )
        )

    def run(self):

        self.ex(
            command=" ".join([
                "run.py",
                "--input",
                self.in_fastq().path,
                "--ref-db",
                self.in_ref_fasta().path,
                "--sample-name",
                "{}.{}".format(self.ref_name, self.sample_name),
                "--output-folder",
                os.path.join(
                    self.base_s3_folder,
                    "align_bwa_{}".format(self.ref_name)
                ),
                "--threads",
                str(self.threads),
                "--temp-folder",
                self.temp_folder
            ])
        )


class FAMLITask(sl.ContainerTask):

    # Inputs: FASTQ and reference database
    in_fastq = None
    in_ref_dmnd = None

    # Parameter: Prefix for output file
    sample_name = sl.Parameter()

    # Parameter: Output folder
    output_folder = sl.Parameter()

    # Parameter: Number of threads for alignment
    threads = sl.Parameter()

    # Parameter: Number of blocks for alignment (each block takes ~6Gb)
    blocks = sl.Parameter(default=5)

    # Parameter: Temporary folder to use on the device
    temp_folder = sl.Parameter()

    # URL of the container
    container = "quay.io/fhcrc-microbiome/famli:v1.1"

    def out_json(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                "{}.json.gz".format(self.sample_name)
            )
        )

    def run(self):

        if self.output_folder[-1] != '/':
            self.output_folder += '/'

        self.ex(
            command=" ".join([
                "famli",
                "align",
                "--input",
                self.in_fastq().path,
                "--sample-name",
                self.sample_name,
                "--ref-db",
                self.in_ref_dmnd().path,
                "--output-folder",
                self.output_folder,
                "--threads",
                str(self.threads),
                "--blocks",
                str(self.blocks),
                "--temp-folder",
                self.temp_folder
            ])
        )
