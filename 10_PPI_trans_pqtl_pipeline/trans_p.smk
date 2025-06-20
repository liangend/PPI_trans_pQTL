configfile: "trans_p_config.yaml"


MODULE=list(range(1, config['Nmodule']+1))
CHRS=list(range(1, config['Nchr']+1))


rule all:
  input:
    config['small_p']

rule z_mat:
  input:
    prot_meta=config['prot_meta'],
  output:
    file_zmat=expand(config['z_DIR']+'mod{{module}}.chr{chr}.txt.gz', chr=CHRS),
    file_sigma=config['sigma_DIR']+'mod{module}.rds'
  params:
    module='{module}',
    dir_sum_stats=config['dir_sum_stats'],
    indep_snp=config['indep_snp'],
    p_thre=config['p_thre']
  script:
    config['tool_DIR']+'make_zmat_interval.R'

rule p:
  input:
    prot_meta=config['prot_meta'],
    file_zmat=config['z_DIR']+'mod{module}.chr{chr}.txt.gz',
    file_sigma=config['sigma_DIR']+'mod{module}.rds'
  output:
    file_p=config['p_DIR']+'p.mod{module}.chr{chr}.rds'
  priority: 10
  params:
    dir_script=config['tool_DIR'], chr='{chr}', module='{module}'
  script:
    config['tool_DIR']+'trans_pco_interval.R'

rule small_p:
  input:
    file_p=expand(config['p_DIR']+'p.mod{module}.chr{chr}.rds', module=MODULE, chr=CHRS)
  output:
    small_p=config['small_p']
  script:
    config['tool_DIR']+'minp.R'

