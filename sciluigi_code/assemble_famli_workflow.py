#!/usr/bin/env python3
"""Assemble a set of FASTQ files, combine the assemblies, and align with FAMLI."""

import os
import re
import uuid
import argparse
import pandas as pd
import sciluigi as sl
from general_tasks import LoadFile
from general_tasks import FastqpTask
from general_tasks import FAMLITask
from assembly_tasks import AssembleMetaSPAdes
from assembly_tasks import AnnotateProkka
from assembly_tasks import IntegrateAssembliesTask
from sra_tasks import ImportSRAFastq


class AssembleFamliWorkflow(sl.WorkflowTask):

    project_name = sl.Parameter()
    input_location = sl.Parameter(default="S3")
    metadata_fp = sl.Parameter()
    base_s3_folder = sl.Parameter()
    sample_column_name = sl.Parameter()
    input_column_name = sl.Parameter()
    metadata_fp_sep = sl.Parameter(default=",")
    assemble_threads = sl.Parameter(default=8)
    assemble_mem = sl.Parameter(default=32000)
    famli_threads = sl.Parameter(default=8)
    famli_mem = sl.Parameter(default=32000)
    metaphlan_threads = sl.Parameter(default=8)
    metaphlan_mem = sl.Parameter(default=32000)
    humann2_threads = sl.Parameter(default=8)
    humann2_mem = sl.Parameter(default=32000)
    aws_job_role_arn = sl.Parameter()
    aws_s3_scratch_loc = sl.Parameter()
    aws_batch_job_queue = sl.Parameter(default="optimal")
    engine = sl.Parameter(default="aws_batch")
    temp_folder = "/scratch"

    def workflow(self):

        # Make sure the project name is alphanumeric
        assert all([
            s.isalnum() or s == "_"
            for s in self.project_name
        ]), "Project name must be alphanumeric"

        # Data can come from either SRA or S3
        assert self.input_location in ["SRA", "S3"]

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
        tasks_metaspades = {}
        tasks_prokka = {}
        tasks_famli = {}

        # Iterate over all of the rows of samples
        for _, r in metadata.iterrows():

            # Get the sample name and the file location
            sample_name = r[self.sample_column_name]
            input_path = r[self.input_column_name]

            # Make a UUID to isolate temp files for this task from any others
            task_uuid = str(uuid.uuid4())[:8]

            # 1. LOAD THE INPUT FILES
            
            if self.input_location == "S3":
                tasks_load_inputs[sample_name] = self.new_task(
                    "load_from_s3_{}".format(sample_name),
                    LoadFile,
                    path=input_path
                )
            elif self.input_location == "SRA":
                assert input_path.startswith("SRR"), input_path

                tasks_load_inputs[sample_name] = self.new_task(
                    "download_from_SRA_{}".format(sample_name),
                    ImportSRAFastq,
                    sra_accession=input_path,
                    base_s3_folder=self.base_s3_folder,
                    input_mount_point="/scratch/{}_get_sra/input/".format(
                        task_uuid),
                    output_mount_point="/scratch/{}_get_sra/output/".format(
                        task_uuid),
                    containerinfo=sl.ContainerInfo(
                        vcpu=1,
                        mem=32000,
                        engine=self.engine,
                        aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                        aws_jobRoleArn=self.aws_job_role_arn,
                        aws_batch_job_queue=self.aws_batch_job_queue,
                        aws_batch_job_prefix=re.sub(
                            '[^a-zA-Z0-9-_]', '_',
                            "get_sra_{}".format(sample_name)
                        ),
                        mounts={
                            "/docker_scratch": {
                                "bind": self.temp_folder,
                                "mode": "rw"
                            }
                        }
                    )
                )
            else:
                raise Exception("Data must be from S3 or SRA")
            
            # 2. CALCULATE FASTQ QUALITY METRICS
            tasks_fastqp[sample_name] = self.new_task(
                "fastqp_{}".format(sample_name),
                FastqpTask,
                summary_path=os.path.join(
                    self.base_s3_folder,
                    "fastqp",
                    sample_name + ".fastqp.tsv"
                ),
                input_mount_point="/scratch/{}_fastqp/input/".format(task_uuid),
                output_mount_point="/scratch/{}_fastqp/output/".format(task_uuid),
                containerinfo=sl.ContainerInfo(
                    vcpu=1,
                    mem=32000,
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

            # 3. ASSEMBLE WITH METASPADES
            tasks_metaspades[sample_name] = self.new_task(
                "metaspades_{}".format(sample_name),
                AssembleMetaSPAdes,
                sample_name=sample_name,
                output_folder=os.path.join(
                    self.base_s3_folder,
                    "metaspades"
                ),
                threads=self.assemble_threads,
                max_mem=int(int(self.assemble_mem)/1000),
                temp_folder=self.temp_folder,
                containerinfo=sl.ContainerInfo(
                    vcpu=int(self.assemble_threads),
                    mem=int(self.assemble_mem),
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    aws_batch_job_prefix=re.sub(
                        '[^a-zA-Z0-9-_]', '_',
                        "metaspades_{}".format(sample_name)
                    ),
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        }
                    }
                )
            )

            # 4. ANNOTATE ASSEMBLIES WITH PROKKA
            tasks_prokka[sample_name] = self.new_task(
                "prokka_{}".format(sample_name),
                AnnotateProkka,
                sample_name=sample_name,
                output_folder=os.path.join(
                    self.base_s3_folder,
                    "prokka"
                ),
                threads=self.assemble_threads,
                temp_folder=self.temp_folder,
                containerinfo=sl.ContainerInfo(
                    vcpu=int(self.assemble_threads),
                    mem=int(self.assemble_mem),
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    aws_batch_job_prefix=re.sub(
                        '[^a-zA-Z0-9-_]', '_',
                        "prokka_{}".format(sample_name)
                    ),
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        }
                    }
                )
            )



        # Assign the output from tasks_load_inputs to the input to tasks_fastqp
        for sample_name in tasks_load_inputs:
            assert sample_name in tasks_fastqp

            if self.input_location == "S3":
                tasks_fastqp[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file
            elif self.input_location == "SRA":
                tasks_fastqp[sample_name].in_fastq = tasks_load_inputs[sample_name].out_fastq

            assert sample_name in tasks_metaspades
            if self.input_location == "S3":
                tasks_metaspades[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file
            elif self.input_location == "SRA":
                tasks_metaspades[sample_name].in_fastq = tasks_load_inputs[sample_name].out_fastq

            assert sample_name in tasks_prokka
            tasks_prokka[sample_name].in_fasta = tasks_metaspades[sample_name].out_fasta

        # 5. COMBINE ASSEMBLIES
        task_integrate_assemblies = self.new_task(
            "integrate_assemblies-{}".format(self.project_name),
            IntegrateAssembliesTask,
            output_prefix=self.project_name,
            output_folder=os.path.join(
                self.base_s3_folder,
                "integrated_assembly"
            ),
            gff_folder=os.path.join(
                self.base_s3_folder,
                "prokka"
            ),
            fastp_folder=os.path.join(
                self.base_s3_folder,
                "prokka"
            ),
            temp_folder=self.temp_folder,
            containerinfo=sl.ContainerInfo(
                vcpu=8,
                mem=120000,
                engine=self.engine,
                aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                aws_jobRoleArn=self.aws_job_role_arn,
                aws_batch_job_queue=self.aws_batch_job_queue,
                aws_batch_job_prefix="integrate_assemblies_{}".format(self.project_name),
                mounts={
                    "/docker_scratch": {
                        "bind": self.temp_folder,
                        "mode": "rw"
                    }
                }
            )
        )

        task_integrate_assemblies.in_fastp_list = [
            t.out_faa for t in tasks_prokka.values()
        ]
        task_integrate_assemblies.in_gff_list = [
            t.out_gff for t in tasks_prokka.values()
        ]

        # 6. ALIGN AGAINST THE ASSEMBLY USING FAMLI
        tasks_famli = {}
        # Iterate over all of the rows of samples
        for _, r in metadata.iterrows():

            # Get the sample name and the file location
            sample_name = r[self.sample_column_name]
            input_path = r[self.input_column_name]

            tasks_famli[sample_name] = self.new_task(
                "famli_{}".format(sample_name),
                FAMLITask,
                sample_name=sample_name,
                output_folder=os.path.join(
                    self.base_s3_folder,
                    "famli"
                ),
                threads=self.famli_threads,
                temp_folder=self.temp_folder,
                containerinfo=sl.ContainerInfo(
                    vcpu=int(self.famli_threads),
                    mem=int(self.famli_mem),
                    engine=self.engine,
                    aws_s3_scratch_loc=self.aws_s3_scratch_loc,
                    aws_jobRoleArn=self.aws_job_role_arn,
                    aws_batch_job_queue=self.aws_batch_job_queue,
                    aws_batch_job_prefix="famli_{}".format(sample_name),
                    mounts={
                        "/docker_scratch": {
                            "bind": self.temp_folder,
                            "mode": "rw"
                        }
                    }
                )
            )
            # Connect the raw FASTQ input
            if self.input_location == "S3":
                tasks_famli[sample_name].in_fastq = tasks_load_inputs[sample_name].out_file
            elif self.input_location == "SRA":
                tasks_famli[sample_name].in_fastq = tasks_load_inputs[sample_name].out_fastq

            # Connect the reference database
            tasks_famli[sample_name].in_ref_dmnd = task_integrate_assemblies.out_daa

        return tasks_famli, tasks_fastqp


if __name__ == "__main__":
    """Assemble a set of FASTQ files, combine the assemblies, and align with FAMLI."""

    parser = argparse.ArgumentParser(description="""
        Assemble a set of FASTQ files, combine the assemblies, and align with FAMLI.""")

    parser.add_argument(
        "--project-name",
        help="Name for entire project (alphanumeric only)",
        required=True
    )

    parser.add_argument(
        "--metadata-fp",
        help = "Table containing metadata",
        required = True
    )

    parser.add_argument(
        "--base-s3-folder",
        help = "Folder to store data in S3",
        required = True
    )

    parser.add_argument(
        "--input-location",
        help = "Location of input data (SRA or S3)",
        default = "S3"
    )

    parser.add_argument(
        "--sample-column-name",
        help = "Column containing the sample name",
        type = str,
        required = True
    )
    parser.add_argument(
        "--metadata-fp-sep",
        help = "Field separator (usually , or <tab>)",
        type = str,
        default = ","
    )
    parser.add_argument(
        "--input-column-name",
        help = "Column containing input location",
        type = str,
        required = True
    )

    parser.add_argument(
        "--assemble-threads",
        type = int,
        default = 4,
        help = "Number of CPUs to use for assembly with metaSPAdes"
    )

    parser.add_argument(
        "--assemble-mem",
        type = int,
        default = 10000,
        help = "Memory to use for assembly with metaSPAdes (MBs)"
    )

    parser.add_argument(
        "--famli-threads",
        type = int,
        default = 4,
        help = "Number of CPUs to use for alignment with FAMLI"
    )

    parser.add_argument(
        "--famli-mem",
        type = int,
        default = 10000,
        help = "Memory to use for alignment with FAMLI (MBs)"
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
        main_task_cls = AssembleFamliWorkflow,
        cmdline_args = [
            "--{}={}".format(
                k.replace("_", "-"),
                v)
            for k, v in args.__dict__.items()
        ]
    )
