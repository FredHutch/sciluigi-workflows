import os
import sciluigi as sl


class HUMAnN2Task(sl.ContainerTask):
    # Input FASTQ file
    in_fastq = None

    # Reference database
    ref_db = sl.Parameter(default="")
    
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
    container = "quay.io/fhcrc-microbiome/humann2:v0.11.1--7"

    def out_json(self):
        # Output is an S3 object
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.output_folder,
                self.sample_name + ".json.gz"
            )
        )

    def run(self):

        if self.output_folder.endswith("/") is False:
            self.output_folder = self.output_folder + "/"

        self.ex(
            command=" ".join([
                "run.py",
                "--input",
                self.in_fastq().path,
                "--sample-name",
                self.sample_name,
                "--output-folder",
                self.output_folder,
                "--ref-db",
                self.ref_db,
                "--threads",
                str(int(self.threads)),
                "--temp-folder",
                self.temp_folder
            ])
        )
