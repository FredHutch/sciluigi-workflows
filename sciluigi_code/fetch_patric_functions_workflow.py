#!/usr/bin/env python3
"""Get all of the functional annotations from PATRIC genomes, link to 16S records."""

import os
import re
import io
import uuid
import boto3
import argparse
import pandas as pd
import sciluigi as sl
from urllib.request import urlopen
from Bio.SeqIO.FastaIO import SimpleFastaParser


def read_tsv_from_s3_as_dataframe(bucket_name, key_name, sep="\t"):
    s3 = boto3.client('s3')
    retr = s3.get_object(Bucket=bucket_name, Key=key_name)

    bytestream = io.BytesIO(retr['Body'].read())
    return pd.read_table(bytestream, sep=sep)


def read_fasta_from_s3(bucket_name, key_name, sep="\t"):
    s3 = boto3.client('s3')
    retr = s3.get_object(Bucket=bucket_name, Key=key_name)

    for header, seq in SimpleFastaParser(io.StringIO(retr['Body'].read().decode('utf-8'))):
        yield header, seq


class TransferFTPtoS3(sl.ContainerTask):
    """Transfer a file from an FTP server to an AWS S3 bucket."""
    # FTP path
    ftp_url = sl.Parameter()
    # S3 path
    s3_url = sl.Parameter()

    # Container with wget
    container = "quay.io/fhcrc-microbiome/python:python-v0.1"

    def out_file(self):
        # File on S3
        return sl.ContainerTargetInfo(
            self,
            self.s3_url
        )

    def run(self):

        output_targets = {
            "out_file": self.out_file()
        }

        self.ex(
            command="wget -O $out_file {}".format(self.ftp_url),
            output_targets=output_targets
        )


class Extract16S(sl.ContainerTask):
    """Extract all of the 16S transcripts from a list of S3 FASTA files."""
    # Input files
    in_fastas = None
    # Folder with all of the transcripts in subfolders
    s3_parent_folder = sl.Parameter()
    # Single flat file with all of the transcripts
    s3_url = sl.Parameter()
    # Temporary folder to use for downloading data inside the Docker container
    temp_folder = sl.Parameter()

    # Container with wget
    container = "quay.io/fhcrc-microbiome/python:python-v0.1"

    def out_file(self):
        # File on S3
        return sl.ContainerTargetInfo(
            self,
            self.s3_url
        )

    def run(self):

        # Save all of the transcripts with " 16S " or " SSU " in the header
        output = {}

        for genome_transcripts in self.in_fastas:

            # Get the transcript names from the FASTA
            bucket, transcript_key = genome_transcripts(
            ).path[5:].split("/", 1)
            
            for header, seq in read_fasta_from_s3(bucket, transcript_key):
                if " 16S " in header or " SSU " in header:
                    header = header.split(" ", 1)[0]
                    assert header not in output, "Duplicated transcript ID, stopping ({})".format(header)

                    output[header] = seq

        # Now write it to S3
        output_bucket, output_key = self.out_file().path[5:].split("/", 1)
        fasta_buffer = io.StringIO()
        for header, seq in output.items():
            fasta_buffer.write(">{}\n{}\n".format(header, seq))

        s3 = boto3.resource('s3')
        s3.Object(output_bucket, output_key).put(Body=fasta_buffer.getvalue())


class ExtractAnnotations(sl.Task):
    """Extract all of annotations and make a TSV keyed by the 16S accession names."""
    # Input files
    in_fastas = None
    in_annotations = None

    # Single flat file with all of the annotations
    s3_url = sl.Parameter()

    # Temporary folder to use for downloading data inside the Docker container
    temp_folder = sl.Parameter()

    def out_file(self):
        # File on S3
        return sl.ContainerTargetInfo(
            self,
            self.s3_url
        )

    def run(self):

        # The final output is going to be a DataFrame with rows as 16S accessions and columns as annotations, values are copy numbers
        output = {}

        for genome_id, genome_transcripts in self.in_fastas.items():
            # Make sure that we also have an annotation for this genome
            assert genome_id in self.in_annotations

            # Get the transcript names from the FASTA
            bucket, transcript_key = genome_transcripts().path[5:].split("/", 1)
            transcript_ids = [
                header.split(" ", 1)[0]
                for header, seq in read_fasta_from_s3(bucket, transcript_key)
                if " 16S " in header or " SSU " in header
            ]

            # Get the annotations for this genome
            bucket, annotation_key = self.in_annotations[genome_id]().path[5:].split("/", 1)
            genome_annotations = read_tsv_from_s3_as_dataframe(bucket, annotation_key)

            # Make a copy number vector
            functional_copy_numbers = genome_annotations["product"].value_counts().to_dict()
            
            # Add the annotations to the list
            for transcript in transcript_ids:
                assert transcript not in output, "Transcript found twice, stopping ({})".format(transcript)
                output[transcript] = functional_copy_numbers

        # Make a single DataFrame
        output = pd.DataFrame(output).fillna(0).T

        # Now write it to S3
        output_bucket, output_key = self.out_file().path[5:].split("/", 1)
        tsv_buffer = io.StringIO()
        output.to_csv(tsv_buffer, sep='\t')
        s3 = boto3.resource('s3')
        s3.Object(output_bucket, output_key).put(Body=tsv_buffer.getvalue())


class FetchPatricFunctions(sl.WorkflowTask):

    s3_folder = sl.Parameter()
    aws_job_role_arn = sl.Parameter()
    aws_s3_scratch_loc = sl.Parameter()
    aws_batch_job_queue = sl.Parameter(default="optimal")
    engine = sl.Parameter(default="aws_batch")
    temp_folder = "/scratch"

    def workflow(self):

        # Make sure that the S3 folder is formatted with the proper prefix
        assert self.s3_folder.startswith("s3://")

        # Parse the bucket and key for the s3 folder for all results
        s3_bucket, s3_prefix = self.s3_folder[5:].split("/", 1)

        # Connect to S3
        s3 = boto3.resource('s3')

        # 1. Get the summary of all genomes
        genome_metadata_fp = os.path.join(s3_prefix, "patric_genome_metadata.tsv")
        
        print("Writing PATRIC genome metadata to s3://{}/{}".format(
            s3_bucket,
            genome_metadata_fp
        ))

        with urlopen("ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata") as fi:
            s3.Bucket(
                s3_bucket
            ).put_object(
                Key=genome_metadata_fp, 
                Body=fi.read()
            )

        # Now read in all of that information as a table
        genome_metadata = read_tsv_from_s3_as_dataframe(s3_bucket, genome_metadata_fp, sep="\t")

        # 2. Fetch the transcripts and annotation files for every genome
        fetch_transcripts_tasks = {}
        fetch_annotation_tasks = {}

        for genome_accession in map(str, genome_metadata.index.values):
            
            fetch_annotation_tasks[genome_accession] = [
                self.new_task(
                    "fetch_patric_annotations_{}".format(genome_accession),
                    TransferFTPtoS3,
                    ftp_url="ftp://ftp.patricbrc.org/genomes/{}/{}.{}".format(
                        genome_accession,
                        genome_accession,
                        suffix
                    ),
                    s3_url=os.path.join(
                        self.s3_folder,
                        genome_accession,
                        "annotation.tsv"
                    ),
                    containerinfo=sl.ContainerInfo(
                        vcpu=1,
                        mem=1000,
                        engine=self.engine,
                        aws_batch_job_poll_sec=120,
                        aws_jobRoleArn=self.aws_job_role_arn,
                        aws_batch_job_queue=self.aws_batch_job_queue,
                        aws_batch_job_prefix=re.sub(
                            '[^a-zA-Z0-9-_]', '_',
                            "fetch_patric_annotations_{}".format(genome_accession)
                        )
                    )
                )
                for suffix in ["PATRIC.pathway.tab", "RefSeq.pathway.tab"]
            ]

            fetch_transcripts_tasks[genome_accession] = [
                self.new_task(
                    "fetch_patric_transcripts_{}".format(genome_accession),
                    TransferFTPtoS3,
                    ftp_url="ftp://ftp.patricbrc.org/genomes/{}/{}.{}".format(
                        genome_accession,
                        genome_accession,
                        suffix
                    ),
                    s3_url=os.path.join(
                        self.s3_folder,
                        genome_accession,
                        "transcripts.frn"
                    ),
                    containerinfo=sl.ContainerInfo(
                        vcpu=1,
                        mem=1000,
                        engine=self.engine,
                        aws_batch_job_poll_sec=120,
                        aws_jobRoleArn=self.aws_job_role_arn,
                        aws_batch_job_queue=self.aws_batch_job_queue,
                        aws_batch_job_prefix=re.sub(
                            '[^a-zA-Z0-9-_]', '_',
                            "fetch_patric_transcripts_{}".format(
                                genome_accession)
                        )
                    )
                )
                for suffix in ["PATRIC.frn", "RefSeq.frn"]
            ]

        # 3. Make a flat file for the 16S records

        extract_all_16S = self.new_task(
            "extract_all_16S",
            Extract16S,
            s3_parent_folder=self.s3_folder,
            s3_url=os.path.join(self.s3_folder, "transcripts.fasta"),
            temp_folder=self.temp_folder,
            containerinfo=sl.ContainerInfo(
                vcpu=1,
                mem=1000,
                engine=self.engine,
                aws_batch_job_poll_sec=120,
                aws_jobRoleArn=self.aws_job_role_arn,
                aws_batch_job_queue=self.aws_batch_job_queue,
                aws_batch_job_prefix="extract_all_16s"
            )
        )

        extract_all_16S.in_fastas = [
            genome_transcript[0].out_file for genome_transcript in fetch_transcripts_tasks.values()
        ]

        # 4. Make a flat file for the annotations
        extract_all_annotations = self.new_task(
            "extract_all_annotations",
            ExtractAnnotations,
            s3_parent_folder=self.s3_folder,
            s3_url=os.path.join(self.s3_folder, "annotations.tsv"),
            temp_folder=self.temp_folder,
            containerinfo=sl.ContainerInfo(
                vcpu=1,
                mem=1000,
                engine=self.engine,
                aws_batch_job_poll_sec=120,
                aws_jobRoleArn=self.aws_job_role_arn,
                aws_batch_job_queue=self.aws_batch_job_queue,
                aws_batch_job_prefix="extract_all_annotations"
            )
        )

        extract_all_annotations.in_fastas = {
            genome_id: genome_transcript[0].out_file
            for genome_id, genome_transcript in fetch_transcripts_tasks.items()
        }
        extract_all_annotations.in_annotations = {
            genome_id: genome_annotation[0].out_file
            for genome_id, genome_annotation in fetch_annotation_tasks.items()
        }

        return extract_all_16S, extract_all_annotations


if __name__ == "__main__":
    """Get all of the functional annotations from PATRIC genomes, link to 16S records."""

    parser = argparse.ArgumentParser(description="""
        Get all of the functional annotations from PATRIC genomes, link to 16S records.""")

    parser.add_argument(
        "--s3-folder",
        help="Folder to fetch all data to, on AWS S3",
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
    args = parser.parse_args()

    sl.run(
        main_task_cls=FetchPatricFunctions,
        cmdline_args=[
            "--{}={}".format(
                k.replace("_", "-"),
                v)
            for k, v in args.__dict__.items()
        ]
    )
