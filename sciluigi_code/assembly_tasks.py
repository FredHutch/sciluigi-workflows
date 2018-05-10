import os
import uuid
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

        if self.output_folder.endswith("/") is False:
            self.output_folder = self.output_folder + "/"

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

    def out_faa(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".fastp.gz"
            )
        )

    def run(self):

        if self.output_folder.endswith("/") is False:
            self.output_folder = self.output_folder + "/"

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


class CheckM(sl.ContainerTask):
    # Input FASTP file of protein sequences
    in_faa = None

    # Sample name
    sample_name = sl.Parameter()
    # Output folder
    output_folder = sl.Parameter()
    # Number of threads to use
    threads = sl.Parameter(default=4)
    # Scratch directory
    temp_folder = sl.Parameter(default="/scratch")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/checkm:checkm-v1.0.11"

    def out_tarball(self):
        # Output is a tarball with all of the results
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".tar.gz"
            )
        )

    def run(self):

        if self.output_folder.endswith("/") is False:
            self.output_folder = self.output_folder + "/"

        input_targets = {
            "faa": self.in_faa()
        }

        output_targets = {
            "tarball": self.out_tarball()
        }

        temp_dir = os.path.join(self.temp_folder, str(uuid.uuid4())[:8])

        self.ex(
            command="echo Checking to see if temp directory ({}) exists && ".format(temp_dir) +
                    "[[ ! -d {} ]] && ".format(temp_dir) +
                    "echo Making temp directory {} && ".format(temp_dir) +
                    "mkdir {} && ".format(temp_dir) +
                    "echo Making temp directory for input files {}/checkm_input && ".format(temp_dir) +
                    "mkdir {}/checkm_input && ".format(temp_dir) +
                    "echo Making temp directory for output files {}/checkm_output && ".format(temp_dir) +
                    "mkdir {}/checkm_output && ".format(temp_dir) +
                    "echo Moving gene FAA file into input directory && " +
                    "mv $faa " + "{}/checkm_input/ && ".format(temp_dir) +
                    "echo Decompressing input file && " +
                    "gunzip {}/checkm_input/* && ".format(temp_dir) +
                    "echo Running checkm && " +
                    "checkm lineage_wf --genes -x faa {}/checkm_input/ {}/checkm_output/ && ".format(
                        temp_dir, temp_dir
                    ) + 
                    "echo Finished running checkm && " +
                    "ls -lhtr {}/checkm_output/ && ".format(temp_dir) +
                    "echo Collecting results && " +
                    "tar cvf {}/output.tar {}/checkm_output && ".format(temp_dir, temp_dir) +
                    "echo Compressing results && " +
                    "gzip {}/output.tar && ".format(temp_dir) +
                    "echo Copying results out of the container && " +
                    "mv {}/output.tar.gz ".format(temp_dir) + "$tarball",
            input_targets=input_targets,
            output_targets=output_targets
            )

