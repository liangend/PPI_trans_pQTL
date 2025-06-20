import pandas as pd
configfile: "h2_cal_config.yaml"
PROT=pd.read_table(config["prot"]).set_index("target",drop=False)

rule all:
  input:
    expand(config['h2_DIR']+'{prot}_h2.log', prot=PROT.index)

## make .sumstats.gz file as ldsc input
rule sumstats:
  input:
    raw_sum_stats=config['raw_sum_stats_DIR']+'{prot}'
  output:
    sumstats=config['sumstats_DIR']+'{prot}.sumstats.gz'
  params:
    snp_list=config['snp_list'],
    n_sample=config['n_sample']
  script:
    config['tool_DIR']+'make_sumstats_ukb.R'

rule h2:
  input:
    sumstats=config['sumstats_DIR']+'{prot}.sumstats.gz'
  output:
    h2_output=config['h2_DIR']+'{prot}_h2.log'
  params:
    ref_annot=config['ref_annot'],
    weight_ld=config['weight_ld'],
    frqfile=config['frqfile'],
    output_prefix=config['h2_DIR']+'{prot}_h2'
  priority: 10
  shell:
    'bash '+config['tool_DIR']+'ldsc_h2.sh {input.sumstats} {params.ref_annot} {params.weight_ld} {params.frqfile} {params.output_prefix}'

