#!/bin/bash
#SBATCH --job-name=trans_pqtl
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=96:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=ana_%j.out
#SBATCH --error=ana_%j.err

mcl /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/bioplex_ppi.txt \
	--abc -I 5.0 \
	-o /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/mcl_result/bioplex_I50

mcl /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/string_ppi.txt \
        --abc -I 5.0 \
	-o /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/mcl_result/string_I50

mcl /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/hippie_ppi.txt \
        --abc -I 5.0 \
	-o /project/xuanyao/jinghui/pqtl/11_prot_complex/ppi_input/mcl_result/hippie_I50

