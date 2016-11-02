import pytest
from click.testing import CliRunner
from parkour import cli
import md5


def file_checksums_equal(file1, file2):
    with open(file1) as f:
        checksum1 = md5.new(f.read()).digest()
    with open(file2) as f:
        checksum2 = md5.new(f.read()).digest()
    return checksum1==checksum2 


def test_trimmed_output():
	runner = CliRunner()
	result = runner.invoke(cli.main, ['-a', 'fastq/s3_1.fastq.gz', '-b', 'fastq/s3_2.fastq.gz', 'trim'])
	print(result.output)
	assert file_checksums_equal('p.s3_1.trim.fastq.gz', 'correct_output/p.s3_1.trim.fastq.gz')
	assert file_checksums_equal('p.s3_2.trim.fastq.gz', 'correct_output/p.s3_2.trim.fastq.gz')
