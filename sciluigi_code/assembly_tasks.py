import os
import sciluigi as sl


class AssembleMetaSPAdes(sl.ContainerTask):
    # Input FASTQ file
    in_fastq = None

    # Sample name
    sample_name = sl.Parameter()
    # Output folder
    output_folder = sl.Parameter()
    # Number of threads to use
    threads = sl.Parameter(default=4)
    # Maximum amount of memory to use (gigabytes)
    max_mem = sl.Parameter(default=10)
    # Scratch directory
    temp_folder = sl.Parameter(default="/scratch")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/metaspades:v3.11.1--7"

    def out_fasta(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".fasta.gz"
            )
        )

    def run(self):

        self.ex(
            command=" ".join([
                "run_metaspades.py",
                "--input",
                self.in_fastq().path,
                "--sample-name",
                self.sample_name,
                "--output-folder",
                self.output_folder,
                "--threads",
                str(int(self.threads)),
                "--max-mem",
                str(int(self.max_mem)),
                "--temp-folder",
                self.temp_folder
            ])
        )


class AnnotateProkka(sl.ContainerTask):
    # Input FASTA file
    in_fasta = None

    # Sample name
    sample_name = sl.Parameter()
    # Output folder
    output_folder = sl.Parameter()
    # Number of threads to use
    threads = sl.Parameter(default=4)
    # Scratch directory
    temp_folder = sl.Parameter(default="/scratch")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/metaspades:v3.11.1--7"

    def out_gff(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".gff.gz"
            )
        )

    def run(self):

        self.ex(
            command=" ".join([
                "run_prokka.py",
                "--input",
                self.in_fasta().path,
                "--sample-name",
                self.sample_name,
                "--output-folder",
                self.output_folder,
                "--threads",
                str(int(self.threads)),
                "--temp-folder",
                self.temp_folder
            ])
        )
