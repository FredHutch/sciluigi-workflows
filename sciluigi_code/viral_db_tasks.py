import os
import sciluigi as sl


class MapVirusesTask(sl.ContainerTask):

    # Inputs: FASTQ and reference database
    in_fastq = None
    in_ref_db_metadata = None
    in_ref_db_dmnd = None

    # Parameter: AWS S3 folder for this project
    base_s3_folder = sl.Parameter()

    # Parameter: Name for output file(s)
    sample_name = sl.Parameter()

    # Parameter: Number of threads for alignment
    threads = sl.Parameter()

    # Parameter: Temporary folder to use on the devide
    temp_folder = sl.Parameter()

    # URL of the container
    container = "quay.io/fhcrc-microbiome/map_viruses:v0.4"

    def out_json(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.base_s3_folder,
                "map_viruses",
                self.sample_name + ".json.gz"
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
                self.temp_folder
            ])
        )
