import os
import sciluigi as sl


class ImportSRAFastq(sl.ContainerTask):
    # Parameter: SRA accession to download data from
    sra_accession = sl.Parameter()

    # Parameter: AWS S3 folder for this project
    base_s3_folder = sl.Parameter()

    # Scratch directory
    scratch_directory = sl.Parameter(default="/scratch")

    # URL of the container
    container = "quay.io/fhcrc-microbiome/get_sra:v0.2"

    def out_fastq(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.base_s3_folder,
                "reads",
                self.sra_accession + ".fastq.gz"
            )
        )

    def run(self):

        self.ex(
            command=" ".join([
                "get_sra.py",
                "--accession",
                self.sra_accession,
                "--output-path",
                self.out_fastq().path,
                "--temp-folder",
                self.scratch_directory
            ])
        )
