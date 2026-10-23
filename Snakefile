COVS=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100]

rule all:
    input:
        "genomes/a.sig.zip",
        "genomes/b.sig.zip",
        #expand("pure/reads/reads-{cov}.sig.zip", cov=COVS),
        expand('pure/mapping.{genome}.csv', genome=['a', 'b']),
        "pure/reads/mf.txt",
        "pure/reads.isect.sig.zip",
        'pure/reads.isect.abund.csv',
        #expand('mixed_90_10/reads/reads-{cov}.sig.zip', cov=COVS),
        "mixed_90_10/reads.isect.sig.zip",
        "mixed_90_10/reads.isect.abund.csv",
        expand('mixed_90_10/mapping.{genome}.csv', genome=['a', 'b']),

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
    input: "genomes/a.fasta"
    output: "pure/reads/reads-{cov}.fasta"
    shell: """
        scripts/make-reads.py {input} -o {output} -C {wildcards.cov}
    """

rule make_mixed_reads_a_b_90_10:
    input:
        "genomes/a.fasta",
        "genomes/b.fasta",
    output:
        'mixed_90_10/reads/reads-{cov}.fasta',
    shell: """
        scripts/make-mixed-reads.py {input} -o {output} -C {wildcards.cov} \
           --mix 9:1
    """

rule sketch_genome:
    input: "genomes/{genome}.fasta",
    output: "genomes/{genome}.sig.zip",
    shell: """
        sourmash sketch dna -p scaled=100,k=21 --name {wildcards.genome} \
           {input} -o {output}
    """
            
rule sketch_reads:
    input: "{dir}/reads/{filename}-{cov}.fasta"
    output: "{dir}/reads/{filename}-{cov}.sig.zip",
    shell: """
        sourmash sketch dna -p scaled=100,k=21,abund --name R-{wildcards.cov} \
           {input} -o {output}
    """

rule make_pure_list:
    input:
        expand("pure/reads/reads-{cov}.sig.zip", cov=COVS),
    output:
        "pure/reads/mf.txt"
    shell: """
        ls -1 {input} > {output}
    """

rule isect_pure:
    input:
        genomes=expand("genomes/{g}.sig.zip", g=['a', 'b']),
        against='pure/reads/mf.txt'
    output:
        "pure/reads.isect.sig.zip"
    shell: """
         ../2026-pancoding-core/scripts/extract-intersection-sigs.py \
             --query-sigs {input.genomes} --metagenomes {input.against} \
             --scaled 100 -k 21 -O {output}
    """

rule make_mixed_list:
    input:
        expand("mixed_90_10/reads/reads-{cov}.sig.zip", cov=COVS),
    output:
        "mixed_90_10/reads/mf.txt"
    shell: """
        ls -1 {input} > {output}
    """

rule isect_mixed_90_10:
    input:
        genomes=expand("genomes/{g}.sig.zip", g=['a', 'b']),
        against='mixed_90_10/reads/mf.txt'
    output:
        "mixed_90_10/reads.isect.sig.zip"
    shell: """
         ../2026-pancoding-core/scripts/extract-intersection-sigs.py \
             --query-sigs {input.genomes} --metagenomes {input.against} \
             --scaled 100 -k 21 -O {output}
    """

rule kmer_est_abund:
    input:
        "{dir}/reads.isect.sig.zip",
    output:
        kmer='{dir}/reads.isect.kmer_stats.csv',
        abund='{dir}/reads.isect.abund.csv',
    shell: """
        ../2026-pancoding-core/scripts/isect-estimate-species-abundances.py \
            {input} --output-kmer {output.kmer} --output-abund {output.abund} \
            --no-unmangle
    """

# NOTE: reads are always 'a' reads, but being mapped to a or b genomes.
rule map_reads:
    input:
        fa='genomes/{genome}.fasta',
        metag='{dir}/reads/reads-{cov}.fasta',
    output:
        bam='{dir}/mapping/R-{cov}.x.{genome}.bam',
        bai='{dir}/mapping/R-{cov}.x.{genome}.bam.bai',
    threads: 8
    shell: """
        minimap2 -ax sr -t {threads} {input.fa:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule bam_to_depth:
    input:
        bam='{dir}/mapping/R-{cov}.x.{genome}.bam',
        bai='{dir}/mapping/R-{cov}.x.{genome}.bam.bai',
        fa='genomes/{genome}.fasta',
    output:
        "{dir}/mapping/R-{cov}.x.{genome}.depth.txt",
    threads: 8
    shell: """
        samtools depth -aa {input.bam:q} {input.fa:q} > {output:q}
    """

# summarize all the highcov mapping depth files -> CSV
rule summarize_depth:
    input:
        expand("{{dir}}/mapping/R-{cov}.x.{{genome}}.depth.txt",
               cov=COVS)
    output:
        '{dir}/mapping.{genome}.csv',
    shell: """
        ../2026-pancoding-core/scripts/summarize-mapping-depth-to-csv.py {input} -o {output}
    """
