import os
import os.path


def calculate_enrich_stats(all_enrich_values_file, output_stats_file, num_mads):
    import sys
    import pandas as pd

    #enrichment_val_df = pd.read_csv(all_enrich_values_file, sep=' ', header=None)
    #enrichment_val_series = enrichment_val_df.iloc[:,-1]
    enrichment_val_series = pd.read_csv(all_enrich_values_file, sep=' ', header=None, usecols=[6], squeeze=True)
    median = enrichment_val_series.median()
    mad = enrichment_val_series.mad()
    num_items=len(enrichment_val_series)
    
    cutoff_value = round(median + (int(num_mads)*mad), 2)
    
    with open(output_stats_file, "w") as out_file:
        out_file.write('Numkmers:' + ('\t') + str(num_items) + '\n')
        out_file.write('MAD' + ('\t') + str(mad) + '\n') 
        out_file.write('Median' + ('\t') + str(median) + '\n')
        out_file.write('Median + ' + str(num_mads) + " mad:" + '\t' + str(cutoff_value) + '\n')
    return(cutoff_value)

workdir: config["WORKING_DIR"]


rule all:
    input:
        # expand("{samplename}/kmer_dbs/{basename}/{kmer_len}/normkcounts.txt", samplename=config["samples"].keys(), basename=config['IPS'], kmer_len=config["KMER_LENS"]),
        #expand("{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.normkcounts.enratio.txt", samplename=config["samples"].keys(), kmer_len=config["KMER_LENS"], pairing=config[
        #    "pairing_config"].keys())
        #expand("{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.enrichratio.{enrichval}.fa", samplename=config["samples"].keys(), kmer_len=config["KMER_LENS"], pairing=config["pairing_config"].keys(), enrichval=config["ENRICH_VAL"])
        #expand("yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract.fa", dataacc=config["your_fav_seqdata"].keys(), samplename=config["samples"].keys(), kmer_len=config["KMER_LENS"], cival=config["CIVAL"], nummads=config["NUM_MADS"], pairing=config["pairing_config"].keys())
        #expand("yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract_singleline.fa", dataacc=config["your_fav_seqdata"].keys(), samplename=config["samples"].keys(), kmer_len=config["KMER_LENS"], cival=config["CIVAL"], nummads=config["NUM_MADS"], pairing=config["pairing_config"].keys())
        expand("yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract_singleline.kmertable", dataacc=config["your_fav_seqdata"].keys(), samplename=config["samples"].keys(), kmer_len=config["KMER_LENS"], cival=config["CIVAL"], nummads=config["NUM_MADS"], pairing=config["pairing_config"].keys())



rule prepare_input:
    input: lambda wildcards: config['samples'][wildcards.samplename][wildcards.basename][wildcards.readnum]
    output: "{samplename}/raw/{basename}/{readnum}_fq.gz"
    threads: 8
    shell:'''
    ln -s $(readlink -f {input}) {output}
    '''

rule prepare_more_input:
    input: lambda wildcards: config['your_fav_seqdata'][wildcards.dataacc]["fq"]
    output: "yf_seqdat/{dataacc}/raw/{dataacc}.fq.gz"
    shell:'''
    ln -s $(readlink -f {input}) {output}
    '''

rule rid_PCR_dups:
    input: "{samplename}/raw/{basename}/{readnum}_fq.gz"
    output: "{samplename}/deduped/{basename}/{readnum}_fq.gz"
    log: "{samplename}/deduped/{basename}/{readnum}_deduping.log"
    threads: 8
    params:
        mem = config["MEM"]
    shell: """
clumpify.sh t={threads} in={input} out={output} dedupe=t subs=0 reorder=f overwrite=t -Xmx{params.mem} 2>{log}
"""

rule trim_adapters_SE:
    input: "{samplename}/deduped/{basename}/{readnum}_fq.gz"
    output: "{samplename}/trim/{basename}/{readnum}_fq.gz"
    threads: 8
    params:
        trimmomatic_path = config["trimmomatic_path"],
        trimmomatic_options = config["trimmomatic_options"],
        trimfasta = config["trimfasta"],
        mem = config["MEM"]
    shell: '''java -XX:ParallelGCThreads={threads} -Xmx{params.mem} -jar {params.trimmomatic_path} SE -phred33 {input} {output} ILLUMINACLIP:{params.trimfasta}:{params.trimmomatic_options}
    '''

# rule pair_reads:
#     input:
#         r1 = "{samplename}/trim/{basename}/fp.fastq.gz",
#         r2 = "{samplename}/trim/{basename}/rp.fastq.gz"
#     output:
#         a = "{samplename}/pear/{basename}/merged.assembled.fastq.gz",
#         uf = "{samplename}/pear/{basename}/merged.unassembled.forward.fastq.gz",
#         ur = "{samplename}/pear/{basename}/merged.unassembled.reverse.fastq.gz",
#         d = "{samplename}/pear/{basename}/merged.discarded.fastq.gz"
#     threads: 8
#     params:
#         PEAR_CONFIG = config["PEAR_CONFIG"]
#     #shell: '''pear -f {input.r1} -r {input.r2} {params.PEAR_CONFIG} -o "{wildcards.samplename}/pear/{wildcards.basename}/merged" 
#     #'''
#     shell: '''pear -f {input.r1} -r {input.r2} {params.PEAR_CONFIG} -o "{wildcards.samplename}/pear/{wildcards.basename}/merged" -j {threads} && (pigz -f -p {threads} "$(dirname {output.a})"/merged.assembled.fastq ; pigz -f -p {threads} "$(dirname {output.a})"/merged.unassembled.forward.fastq ; pigz -f -p {threads} "$(dirname {output.a})"/merged.unassembled.reverse.fastq; pigz -f -p {threads} "$(dirname {output.a})"/merged.discarded.fastq) 
#     '''



rule feed_kmc:
    input: "{samplename}/trim/{basename}/{readnum}_fq.gz"
    output: "{samplename}/kmer_dbs/{basename}_{readnum}_fq.gz"
    shell:'''
    rm -rf "$(dirname {output})"/tmp
    mkdir -p "$(dirname {output})"/tmp
    ln -s $(readlink -f {input}) {output}
    '''

rule generate_kmers:
    input: "{samplename}/kmer_dbs/{basename}_{readnum}_fq.gz"
    output: "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.kmc_suf", "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.kmc_pre"
    params:
        mem = config["MEM"],
        dbname = "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts"
    threads: 8
    shell: '''
    kmc -k{wildcards.kmer_len} -m{params.mem} -sm -t{threads} -ci{wildcards.cival} -cs100000000 {input} {params.dbname} "$(dirname {input})"/tmp
    '''

rule compute_total_bases:
    input: "{samplename}/trim/{basename}/{readnum}_fq.gz"
    output: 
        readlens="{samplename}/kmer_dbs/{basename}_{readnum}/readlens.csv"
    #params:
    #    fq_list=lambda wildcards: " ".join(["{s}/trim/{b}/{type}.fastq.gz".format(s=wildcards.samplename, b=wildcards.basename,type=t) for t in config["types_to_keep"]])
    shell:'''
    zcat {input} | sed -e 's/\t/ /g' | awk -F '[ ]' 'BEGIN{{OFS=" ";}}(NR%4==1){{x=$0; y=(NR+3)/4; getline; printf("%d\\n",length($1));}}' | awk -F '/n' '{{sum+=$1}} END {{print sum}}' > {output.readlens}
    '''

rule dump_kmers:
    input: "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.kmc_pre", "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.kmc_suf"
    output:
        #dump = "{samplename}/kmer_dbs/{basename}/{kmer_len}/kcounts.dump",
        dump_sort = "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.dump.sort"
    params:
        dbname = "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts"
    shell: '''kmc_tools transform {params.dbname} dump -s {output.dump_sort}
    '''

rule norm_kcounts:
    input:
        readlens="{samplename}/kmer_dbs/{basename}_{readnum}/readlens.csv",
        dump_sort = "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/kcounts.dump.sort"

    output: "{samplename}/kmer_dbs/{basename}_{readnum}/{kmer_len}/{cival}/normkcounts.txt"
    shell:'''
    mylen="$(cat {input.readlens})"
    awk -v x=$mylen '{{print$0, $2/x}}' {input.dump_sort} > {output}
    '''
rule join_enrichment_ratio:
    input: 
        in1 = lambda wildcards: "{samplename}/kmer_dbs/{pairing}/{kmer_len}/{cival}/normkcounts.txt".format(samplename=wildcards.samplename, kmer_len=wildcards.kmer_len, cival=wildcards.cival, pairing=config["pairing_config"][wildcards.pairing][0]),
        in2 = lambda wildcards: "{samplename}/kmer_dbs/{pairing}/{kmer_len}/{cival}/normkcounts.txt".format(samplename=wildcards.samplename, kmer_len=wildcards.kmer_len, cival=wildcards.cival, pairing=config["pairing_config"][wildcards.pairing][1])
    output:
        o1 = "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.normkcounts.enratio.txt"
    shell: '''join --check-order {input.in1} {input.in2} | awk -F '[ ]' 'BEGIN{{OFS=" ";}}{{print $0, $2/$4, $3/$5;}}' | sort -r -n -k 7 > {output.o1}
    '''

# rule chop_enratio:
#     input: "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.normkcounts.enratio.txt"
#     output: "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.enrichratio.{enrichval}.txt"
#     shell:'''
#     awk -v enrich={wildcards.enrichval} '{{if ($7 >= enrich) print$0}}' {input} > {output}
#     '''

rule enrichment_cutoff:
    input: "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.normkcounts.enratio.txt"
    output: 
        o1 = "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.madx{nummads}.txt",
        o2 = "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.enrichmentstats_madx{nummads}.txt"
    threads: 1
    run:
        print('started enrichment_cutoff')
        cutoff_value=calculate_enrich_stats(str(input),str(output.o2),wildcards.nummads)
        #cutoff_value=calculate_enrich_stats(str(input),wildcards.nummads)
        print('finished enrichment_cutoff')
        shell("awk -v enrich={{cutoff_value}}".format(cutoff_value=cutoff_value)+" '{{if ($7 >= enrich) print$0}}' {input} > {output.o1}")
# rule generate_enrichment_statsfile:
#     input: 
#         in1 = "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.normkcounts.enratio.txt"
#     output: 
#         # o1 = "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.stdev_x_{numstd}.txt"
#         o2 = "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.enrichment.cutoffvalues.txt"
#     run:
#         cutoff_value=calculate_enrich_stats(input,{output.o2},wildcards.numstd)
#         # shell("awk -v enrich={cutoff_value}".format(cutoff_value=cutoff_value)+"'{{if ($7 >= enrich) print$0}}' {input} > {output.o1}")
 
# rule chop_enratio:
#     input: 
#         all_in = "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.normkcounts.enratio.txt"
#         cutoffs = "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.enrichment.cutoffvalues.txt"
#     output: "{samplename}/kmer_dbs/joined/{kmer_len}/{pairing}.enrichratio.{std_cutoff}.txt"
#     run:
#     for line in input.cutoffs:
#         if line.startswitch('enrichment_cutoff'):
#             cutoff1=line.strip.split('\t')[0].split(' ')[1]
#             cutoff2=line.strip.split('\t')[1].split(' ')[1]
#             cutoff3=line.strip.split('\t')[2].split(' ')[1]
#     shell:'''
#     awk -v enrich={wildcards.enrichval} '{{if ($7 >= enrich) print$0}}' {input} > {output}
#     '''

rule enrichment_to_fa:
    input: "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.madx{nummads}.txt"

    output: 
        o1="{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.enrichratio.madx{nummads}.fa",
        o2="{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.enrichratio.madx{nummads}_clipname.fa"
    shell: '''
    awk -F '[ ]' 'BEGIN{{OFS="\\n"}}{{n=NR; x=">"n"_"$7; print x, $1}}' {input}>{output.o1}
    awk -F '[ ]' 'BEGIN{{OFS="\\n"}}{{n=NR; x=">"n; print x, $1}}' {input}>{output.o2}
    '''

rule bbduk:
    input:
        in1 = "yf_seqdat/{dataacc}/raw/{dataacc}.fq.gz",
        ref = "{samplename}/kmer_dbs/joined/{kmer_len}/{cival}/{pairing}.enrichratio.madx{nummads}_clipname.fa"
    output:
        m = "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract.fa",
        s = "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract.stats.txt"
    shell:'''bbduk.sh in={input.in1} ref={input.ref} outm={output.m} out=/dev/null stats={output.s} overwrite=t k={wildcards.kmer_len} rcomp=t maskmiddle=f mink=-1 rename=t
    '''

rule singleline_fa:
    input: "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract.fa"
    output: "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract_singleline.fa"
    shell:'''
    awk 'BEGIN{{ORS=""; getline; print $0"\\n"}}(/^>/){{print "\\n"$0"\\n"; next;}}{{print $0}}END{{print "\\n"}}' {input} > {output}
    '''    

rule kmer_table:
    input: "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract_singleline.fa"
    output: "yf_seqdat/{dataacc}/{samplename}/{kmer_len}/{cival}/madx{nummads}/{pairing}_extract_singleline.kmertable"
    shell:'''
    cat {input} | sed -e 's/\t/ /' | awk -F '[ ]' 'BEGIN{{N=0; OFS="\t"}}(NR%2==1){{if(NF>2){{N+=1; readid=$1; kmers=$3; getline; print N, readid, length($0), kmers;}}}}' > {output}
    '''  


##run as: snakemake -s snakefile_jstrain_lr.py --configfile jstrain_lr_config.py -pr --cores <cpus>
