#!/usr/bin/env python3
"""Map a set of samples against a viral reference database."""

import os
import argparse
import pandas as pd
import sciluigi as sl
from general_tasks import LoadFile
from viral_db_tasks import MapVirusesTask
from sra_tasks import ImportSRAFastq


class MapVirusesWorkflow(sl.WorkflowTask):

    # If input files are supplied from SRA, fetch them
    # Align the WGS reads against a viral database

    ref_db_metadata = sl.Parameter()
    ref_db_dmnd = sl.Parameter()
    metadata_fp = sl.Parameter()
    base_s3_folder = sl.Parameter()
    sample_column_name = sl.Parameter()
    metadata_fp_sep = sl.Parameter(default=",")
    input_column_name = sl.Parameter()
    input_location = sl.Parameter()
    align_threads = sl.Parameter(default=4)
    align_mem = sl.Parameter(default=10000)
    aws_job_role_arn = sl.Parameter()
    aws_s3_scratch_loc = sl.Parameter()
    aws_batch_job_queue = sl.Parameter(default="optimal")
    engine = sl.Parameter(default="aws_batch")
    temp_folder = "/scratch"

    def workflow(self):

        # Input files are either located in SRA or AWS S3
        assert self.input_location in ["SRA", "S3"]

        # Read in the metadata sheet
        metadata = pd.read_table(self.metadata_fp, sep=self.metadata_fp_sep)
        for col_name in [self.input_column_name, self.sample_column_name]:
            assert col_name in metadata.columns, "{} not found in {}".format(
                col_name, self.metadata_fp
            )
            # Make sure that all samples and files are unique
            assert metadata[col_name].unique().shape[0] == metadata.shape[0]

        # Make tasks that will make sure the reference databases exist
        ref_db_dmnd = self.new_task(
            "load_ref_db_dmnd",
            LoadFile,
            path=self.ref_db_dmnd
        )
        ref_db_metadata = self.new_task(
            "load_ref_db_metadata",
            LoadFile,
            path=self.ref_db_metadata
        )

        # Keep track of all of the jobs for getting the input files
        tasks_load_inputs = {}

        # Keep track of all of the jobs for aligning against the viral database
        tasks_map_viruses = {}

        # Iterate over all of the rows of samples
        for ix, r in metadata.iterrows():

            # Get the sample name and the file location
            sample_name = r[self.sample_column_name]
            input_path = r[self.input_column_name]

            # If the inputs are on SRA, execute jobs that will download them
            if self.input_location == "SRA":

                tasks_load_inputs[sample_name] = self.new_task(
                    "download_from_sra_{}".format(sample_name),
                    ImportSRAFastq,
                    sra_accession=input_path,
                    base_s3_folder=self.base_s3_folder,
                    containerinfo=sl.ContainerInfo(
                        vcpu=1,
                        mem=4096,
                        engine=self.engine,
                        aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                        aws_jobRoleArn=self.aws_job_role_arn,
                        aws_batch_job_queue=self.aws_batch_job_queue,
                        mounts={
                            "/docker_scratch": {
                                "bind": self.temp_folder,
                                "mode": "rw"
                            }
                        }
                    )
                )
            else:
                # Make sure the file exists on S3
                assert self.input_location == "S3"
                tasks_load_inputs[sample_name] = self.new_task(
                    "load_from_s3_{}".format(sample_name),
                    LoadFile,
                    path=input_path
                )

            # Make a task to align the reads, wherever they came from
            tasks_map_viruses[sample_name] = self.new_task(
                "map_viruses_{}".format(sample_name),
                MapVirusesTask,
                base_s3_folder=self.base_s3_folder,
                sample_name=sample_name,
                threads=self.align_threads,
                temp_folder=self.temp_folder,
                containerinfo=sl.ContainerInfo(
                    vcpu=int(self.align_threads),
                    mem=int(self.align_mem),
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        }
                    }
                )
            )
        # Assign the output from tasks_load_inputs to the input to tasks_map_viruses
        for sample_name in tasks_load_inputs:
            assert sample_name in tasks_map_viruses

            # Assign the input for the reference database
            tasks_map_viruses[sample_name].in_ref_db_dmnd = ref_db_dmnd.out_file
            tasks_map_viruses[sample_name].in_ref_db_metadata = ref_db_metadata.out_file
            tasks_map_viruses[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file

        return tasks_map_viruses


if __name__ == "__main__":
    """Map WGS data against a database of viruses."""

    parser = argparse.ArgumentParser(description="""
        Get FASTQ WGS data and map against a database of viruses.""")

    parser.add_argument(
        "--ref-db-metadata",
        help="Location of viral database metadata file",
        required=True
    )

    parser.add_argument(
        "--ref-db-dmnd",
        help="Location of viral database DMND file",
        required=True
    )

    parser.add_argument(
        "--metadata-fp",
        help="Table containing metadata",
        required=True
    )

    parser.add_argument(
        "--base-s3-folder",
        help="Folder to store data in S3",
        required=True
    )

    parser.add_argument(
        "--sample-column-name",
        help="Column containing the sample name",
        type=str,
        required=True
    )
    parser.add_argument(
        "--metadata-fp-sep",
        help="Field separator (usually , or <tab>)",
        type=str,
        default=","
    )
    parser.add_argument(
        "--input-column-name",
        help="Column containing input location",
        type=str,
        required=True
    )
    parser.add_argument(
        "--input-location",
        type=str,
        required=True,
        help="Specify the location of the inputs [SRA | S3]"
    )

    parser.add_argument(
        "--align-threads",
        type=int,
        default=4,
        help="Number of CPUs to use for alignment"
    )

    parser.add_argument(
        "--align-mem",
        type=int,
        default=10000,
        help="Memory to use for alignment (MBs)"
    )

    parser.add_argument(
        "--aws-job-role-arn",
        help="Job Role ARN to use with AWS Batch",
        type=str
    )

    parser.add_argument(
        "--aws-batch-job-queue",
        help="Job Queue to use with AWS Batch",
        type=str,
        default="optimal"
    )

    parser.add_argument(
        "--aws-s3-scratch-loc",
        help="S3 bucket to use for scratch files",
        type=str
    )

    parser.add_argument(
        "--engine",
        help="Execution engine",
        type=str,
        default="aws_batch"
    )

    parser.add_argument(
        "--workers",
        help="Number of workers to use for parallel execution",
        type=int,
        default=500
    )

    args = parser.parse_args()

    # Either specify the SRA or S3
    assert args.input_location in ["SRA", "S3"]

    assert os.path.exists(args.metadata_fp)

    sl.run(
        main_task_cls=MapVirusesWorkflow,
        cmdline_args=[
            "--{}={}".format(
                k.replace("_", "-"), 
                v)
            for k, v in args.__dict__.items()
        ]
    )
