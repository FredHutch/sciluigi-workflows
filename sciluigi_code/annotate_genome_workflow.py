#!/usr/bin/env python3
"""Annotate a single genome."""

import os
import argparse
import sciluigi as sl
from general_tasks import LoadFile
from assembly_tasks import AnnotateProkka, CheckM


class AnnotateGenomeWorkflow(sl.WorkflowTask):
    """ {genome sequence } -> [ prodigal / prokka ] -> {called peptides} -> [checkM] -> {completeness and taxID}"""

    genome_fasta = sl.Parameter()
    genome_name = sl.Parameter()
    checkm_memory = sl.Parameter(default=64000)
    checkm_threads = sl.Parameter(default=8)
    base_s3_folder = sl.Parameter()
    aws_job_role_arn = sl.Parameter()
    aws_s3_scratch_loc = sl.Parameter()
    aws_batch_job_queue = sl.Parameter(default="optimal")
    engine = sl.Parameter(default="aws_batch")
    temp_folder = "/scratch"

    def workflow(self):

        # Load the input file
        genome_fasta = self.new_task(
            "load_genome_fasta",
            LoadFile,
            path=self.genome_fasta
        )

        # Run Prokka
        annotate_prokka = self.new_task(
            "annotate_prokka_{}".format(self.genome_name),
            AnnotateProkka,
            sample_name=self.genome_name,
            output_folder=os.path.join(self.base_s3_folder, "prokka"),
            threads=self.checkm_threads,
            temp_folder=self.temp_folder,
            containerinfo=sl.ContainerInfo(
                vcpu=int(self.checkm_threads),
                mem=int(self.checkm_memory),
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

        # Link the file for prokka annotation
        annotate_prokka.in_fasta = genome_fasta.out_file

        # Run CheckM
        checkm = self.new_task(
            "checkm_{}".format(self.genome_name),
            CheckM,
            sample_name=self.genome_name,
            output_folder=os.path.join(self.base_s3_folder, "checkm"),
            threads=8,
            temp_folder=self.temp_folder,
            containerinfo=sl.ContainerInfo(
                vcpu=int(8),
                mem=int(64000),
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

        # Link the protein coding sequences from prokka into the inputs for checkm
        checkm.in_faa = annotate_prokka.out_faa

        return checkm


if __name__ == "__main__":
    """Annotate a bacterial genome."""

    parser = argparse.ArgumentParser(description="""
        Annotate a bacterial genome and run checkm.""")

    parser.add_argument(
        "--genome-fasta",
        help="Location of FASTA file",
        required=True
    )

    parser.add_argument(
        "--genome-name",
        help="Name of the genome (prefix for output files)",
        required=True
    )

    parser.add_argument(
        "--base-s3-folder",
        help="Folder to store data in S3",
        required=True
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

    sl.run(
        main_task_cls=AnnotateGenomeWorkflow,
        cmdline_args=[
            "--{}={}".format(
                k.replace("_", "-"),
                v)
            for k, v in args.__dict__.items()
        ]
    )
