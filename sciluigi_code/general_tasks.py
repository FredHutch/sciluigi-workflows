import os
import sciluigi as sl

class LoadFile(sl.ExternalTask):
    path = sl.Parameter()

    def out_file(self):
        return sl.ContainerTargetInfo(self, self.path)


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

    # Parameter: Temporary folder to use on the devide
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
