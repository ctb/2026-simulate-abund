COVS=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100]

rule all:
    input:
        "genomes/a.sig.zip",
        "genomes/b.sig.zip",
        expand("reads/a.reads-{cov}.sig.zip", cov=COVS),
        expand("mapping/R-{cov}.x.a.depth.txt", cov=COVS),
        expand("mapping/R-{cov}.x.b.depth.txt", cov=COVS),
        expand('reads.mapping.{genome}.csv', genome=['a', 'b']),
        "reads/mf.csv",
        "reads/mf.txt",
        "reads.isect.sig.zip",
        'reads.isect.abund.csv',

rule make_a_fasta:
    output: "genomes/a.fasta"
    shell: """
        scripts/make-random-cds.py A -o {output}
    """

rule mutate_a_to_b:
    input: "genomes/a.fasta"
    output: "genomes/b.fasta"
    shell: """
        scripts/mutate-cds.py {input} -m 3 -o {output}
    """
            
rule make_reads_a:
    input: "genomes/{f}.fasta"
    output: "reads/{f}.reads-{cov}.fasta"
    shell: """
        scripts/make-reads.py {input} -o {output} -C {wildcards.cov}
    """

rule sketch_genome:
    input: "genomes/{genome}.fasta"
    output: "genomes/{genome}.sig.zip",
    shell: """
        sourmash sketch dna -p scaled=100,k=21 --name {wildcards.genome} \
           {input} -o {output}
    """
            
rule sketch_reads:
    input: "reads/{f}.reads-{cov}.fasta"
    output: "reads/{f}.reads-{cov}.sig.zip",
    shell: """
        sourmash sketch dna -p scaled=100,k=21,abund --name R-{wildcards.cov} \
           {input} -o {output}
    """

rule mf:
    input:
        expand("reads/a.reads-{cov}.sig.zip", cov=COVS),
    output:
        "reads/mf.csv"
    shell: """
        sourmash sig collect {input} -o {output} -F csv --abspath
    """

rule make_list:
    input:
        expand("reads/a.reads-{cov}.sig.zip", cov=COVS),
    output:
        "reads/mf.txt"
    shell: """
        ls -1 {input} > {output}
    """

rule isect:
    input:
        genomes=expand("genomes/{g}.sig.zip", g=['a', 'b']),
        against='reads/mf.txt'
    output:
        "reads.isect.sig.zip"
    shell: """
         ../2026-pancoding-core/scripts/extract-intersection-sigs.py \
             --query-sigs {input.genomes} --metagenomes {input.against} \
             --scaled 100 -k 21 -O {output}
    """

rule kmer_est_abund:
    input:
        "reads.isect.sig.zip"
    output:
        kmer='reads.isect.kmer_stats.csv',
        abund='reads.isect.abund.csv',
    shell: """
        ../2026-pancoding-core/scripts/isect-estimate-species-abundances.py \
            {input} --output-kmer {output.kmer} --output-abund {output.abund} \
            --no-unmangle
    """

# index
rule index_genomes:
    input:
        "genomes/{genome}.fasta",
    output:
        "genomes/{genome}.1.bt2",
    shell: """
        bowtie2-build {input} genomes/{wildcards.genome}
    """

# NOTE: reads are always 'a' reads, but being mapped to a or b genomes.
rule map_reads:
    input:
        fa='genomes/{genome}.fasta',
        metag='reads/a.reads-{cov}.fasta',
        idx="genomes/{genome}.1.bt2",
    output:
        sam=temporary("mapping/R-{cov}.x.{genome}.sam"),
        bam='mapping/R-{cov}.x.{genome}.bam',
        bai='mapping/R-{cov}.x.{genome}.bam.bai',
    threads: 8
    shell: """
        bowtie2 -f -p {threads} --sensitive \
            -x genomes/{wildcards.genome} \
            -U {input.metag} -S {output.sam}
        samtools view -b -F 4 {output.sam} | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule bam_to_depth:
    input:
        bam='mapping/R-{cov}.x.{genome}.bam',
        bai='mapping/R-{cov}.x.{genome}.bam.bai',
        fa='genomes/{genome}.fasta',
    output:
        "mapping/R-{cov}.x.{genome}.depth.txt",
    threads: 8
    shell: """
        samtools depth -aa {input.bam:q} {input.fa:q} > {output:q}
    """

# summarize all the highcov mapping depth files -> CSV
rule summarize_depth:
    input:
        expand("mapping/R-{cov}.x.{{genome}}.depth.txt",
               cov=COVS)
    output:
        'reads.mapping.{genome}.csv',
    shell: """
        ../2026-pancoding-core/scripts/summarize-mapping-depth-to-csv.py {input} -o {output} --exclude-ends 1
    """
