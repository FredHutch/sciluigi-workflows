import os
import sciluigi as sl


class MapVirusesTask(sl.ContainerTask):

    # Inputs: FASTQ and reference database
    in_fastq = None
    in_ref_db_metadata = None
    in_ref_db_dmnd = None

    # Parameter: AWS S3 folder for output files
    output_folder = sl.Parameter()

    # Parameter: Name for output file(s)
    sample_name = sl.Parameter()

    # Parameter: Number of threads for alignment
    threads = sl.Parameter()

    # Parameter: Temporary folder to use on the devide
    temp_folder = sl.Parameter()

    # URL of the container
    container = "quay.io/fhcrc-microbiome/map_viruses:v0.7"

    def out_json(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".json.gz"
            )
        )

    def out_sam(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".sam.gz"
            )
        )

    def run(self):

        self.ex(
            command=" ".join([
                "map_viruses.py",
                "--input",
                self.in_fastq().path,
                "--metadata",
                self.in_ref_db_metadata().path,
                "--ref-db",
                self.in_ref_db_dmnd().path,
                "--output-path",
                self.out_json().path,
                "--threads",
                str(self.threads),
                "--temp-folder",
                self.temp_folder,
                "--keep-alignments",
            ])
        )


class VirFinderTask(sl.ContainerTask):
    """Run the VirFinder tool on a set of contigs."""

    # Inputs: FASTA
    in_fasta = None

    # Parameter: AWS S3 folder for this project
    base_s3_folder = sl.Parameter()

    # Parameter: Name for output file(s)
    sample_name = sl.Parameter()

    # Use this to specify a different mount point for temporary files, if needed
    input_mount_point = sl.Parameter(default="/mnt/input/")
    output_mount_point = sl.Parameter(default="/mnt/output/")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/virfinder:v1.1--0"

    def out_tsv(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.base_s3_folder,
                "virfinder",
                self.sample_name + ".tsv"
            )
        )


    def run(self):

        input_targets={
            "input_fasta": self.in_fasta()
        }

        output_targets={
            "output_tsv": self.out_tsv()
        }

        self.ex(
            command="run_virfinder.Rscript $input_fasta $output_tsv",
            input_targets= input_targets,
            output_targets= output_targets,
            input_mount_point= self.input_mount_point,
            output_mount_point= self.output_mount_point,
        )
