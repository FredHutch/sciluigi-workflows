#!/usr/bin/env python3
"""Analyze a set of FASTQ files with HUMAnN2."""

import os
import re
import uuid
import argparse
import pandas as pd
import sciluigi as sl
from general_tasks import LoadFile
from general_tasks import FastqpTask
from biobakery_tasks import HUMAnN2Task


class HUMAnN2Workflow(sl.WorkflowTask):

    metadata_fp = sl.Parameter()
    base_s3_folder = sl.Parameter()
    sample_column_name = sl.Parameter()
    input_column_name = sl.Parameter()
    metadata_fp_sep = sl.Parameter(default=",")
    humann2_threads = sl.Parameter(default=8)
    humann2_mem = sl.Parameter(default=32000)
    humann2_ref_db = sl.Parameter()
    aws_job_role_arn = sl.Parameter()
    aws_s3_scratch_loc = sl.Parameter()
    aws_batch_job_queue = sl.Parameter(default="optimal")
    engine = sl.Parameter(default="aws_batch")
    temp_folder = "/scratch"

    def workflow(self):

        # Read in the metadata sheet
        metadata = pd.read_table(self.metadata_fp, sep=self.metadata_fp_sep)

        for col_name in [self.input_column_name, self.sample_column_name]:
            assert col_name in metadata.columns, "{} not found in {}".format(
                col_name, self.metadata_fp
            )
            # Make sure that all samples and files are unique
            assert metadata[col_name].unique().shape[0] == metadata.shape[0]

        # Keep track of the jobs for each step, for each sample
        tasks_load_inputs = {}
        tasks_fastqp = {}
        tasks_humann = {}

        # Iterate over all of the rows of samples
        for ix, r in metadata.iterrows():

            # Get the sample name and the file location
            sample_name = r[self.sample_column_name]
            input_path = r[self.input_column_name]

            # Make a UUID to isolate temp files for this task from any others
            task_uuid = str(uuid.uuid4())[:8]

            # 1. LOAD THE INPUT FILES

            tasks_load_inputs[sample_name] = self.new_task(
                "load_from_s3_{}".format(sample_name),
                LoadFile,
                path=input_path
            )

            # 2. CALCULATE FASTQ QUALITY METRICS
            tasks_fastqp[sample_name] = self.new_task(
                "fastqp_{}".format(sample_name),
                FastqpTask,
                summary_path=os.path.join(
                    self.base_s3_folder,
                    "fastqp",
                    sample_name + ".fastqp.tsv"
                ),
                input_mount_point="/scratch/{}_fastqp/input/".format(
                    task_uuid),
                output_mount_point="/scratch/{}_fastqp/output/".format(
                    task_uuid),
                containerinfo=sl.ContainerInfo(
                    vcpu=1,
                    mem=10000,
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    aws_batch_job_prefix=re.sub(
                        '[^a-zA-Z0-9-_]', '_',
                        "fastqp_{}".format(sample_name)
                    ),
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        }
                    }
                )
            )

            # 3. ANALYZE WITH HUMANn2
            tasks_humann[sample_name] = self.new_task(
                "humann2_{}".format(sample_name),
                HUMAnN2Task,
                sample_name=sample_name,
                output_folder=os.path.join(
                    self.base_s3_folder,
                    "humann2"
                ),
                threads=self.humann2_threads,
                ref_db=self.humann2_ref_db,
                temp_folder=self.temp_folder,
                containerinfo=sl.ContainerInfo(
                    vcpu=int(self.humann2_threads),
                    mem=int(self.humann2_mem),
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    aws_batch_job_prefix=re.sub(
                        '[^a-zA-Z0-9-_]', '_',
                        "humann2_{}".format(sample_name)
                    ),
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        },
                        "/refdbs": {
                            "bind": "/refdbs",
                            "mode": "ro"
                        }
                    }
                )
            )

        # Assign the output from tasks_load_inputs to the input to tasks_fastqp
        for sample_name in tasks_load_inputs:
            assert sample_name in tasks_fastqp
            tasks_fastqp[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file

            assert sample_name in tasks_humann
            tasks_humann[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file

        return tasks_fastqp, tasks_humann


if __name__ == "__main__":
    """Assemble a set of FASTQ files, combine the assemblies, and align with FAMLI."""

    parser = argparse.ArgumentParser(description="""
        Assemble a set of FASTQ files, combine the assemblies, and align with FAMLI.""")

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
        "--humann2-threads",
        type=int,
        default=4,
        help="Number of CPUs to use for assembly with metaSPAdes"
    )

    parser.add_argument(
        "--humann2-mem",
        type=int,
        default=10000,
        help="Memory to use for assembly with metaSPAdes (MBs)"
    )

    parser.add_argument(
        "--humann2-ref-db",
        type = str,
        help = "Reference database for HUMAnN2"
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

    assert os.path.exists(args.metadata_fp)

    sl.run(
        main_task_cls = HUMAnN2Workflow,
        cmdline_args=[
            "--{}={}".format(
                k.replace("_", "-"),
                v)
            for k, v in args.__dict__.items()
        ]
    )
