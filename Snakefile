import csv
import sys
import hashlib
from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

with open("data/10x.csv") as f_in:
    reader = csv.DictReader(f_in)
    DATA = list(reader)
    URLS = {Path(row["URL"]).name: row["URL"] for row in DATA}
    FILES = [str(Path("data") / row["Group"] / Path(row["URL"]).name) for row in DATA]

rule all:
    input: "analysis/vdj_v1_hs_pbmc3_b_26x91/vdj_v1_hs_pbmc3_b_26x91.mri.tgz"

rule get_data_all:
    input: FILES

rule get_file:
    output: "data/{group}/{file}"
    input: lambda w: HTTP.remote(URLS[w.file], keep_local=True)
    params: url=lambda w: URLS[w.file]
    run:
        md5_exp = [row for row in DATA if row["URL"] == params.url][0]["md5sum"]
        md5_obs = md5(input[0])
        if not md5_obs == md5_exp:
            raise UserWarning(
                "MD5 mismatch on %s: expected %s, received %s" % (fname, md5_exp, md5_obs))
        shell("mv {input} {output}")

rule reference:
    output: "data/Reference/{name}/reference.json"
    input: "data/Reference/{name}.tar.gz"
    shell: "tar xzf {input} -C data/Reference"

rule cr_vdj:
    output: "analysis/vdj_v1_hs_pbmc3_b_26x91/vdj_v1_hs_pbmc3_b_26x91.mri.tgz"
    input:
        fastqdir="analysis/vdj_v1_hs_pbmc3_b_26x91_fastqs",
        reference="data/Reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0/reference.json"
    params:
        refdir="data/Reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"
    threads: 20
    shell:
        """
            cd analysis
            cellranger vdj \
                --fastqs=../{input.fastqdir} \
                --localcores={threads} \
                --id=vdj_v1_hs_pbmc3_b_26x91 \
                --reference=../{params.refdir} \
                --description='Human PBMC - Ig enrichment from amplified cDNA (v1.0, 26x91)'
        """

rule analysis_fastqs:
    output: directory("analysis/vdj_v1_hs_pbmc3_b_26x91_fastqs")
    input: "data/Input/vdj_v1_hs_pbmc3_b_fastqs.tar"
    shell: "tar xf {input} -C analysis"

rule cr_sitecheck:
    output: "analysis/sitecheck.txt"
    shell: "cellranger sitecheck > {output}"

rule cr_testrun:
    output: "tiny/tiny.mri.tgz"
    shell: "cellranger testrun --id=tiny"

# https://stackoverflow.com/a/3431838
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
