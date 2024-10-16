#!/bin/bash
#SBATCH --job-name=trans_eqtl
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10GB
#SBATCH --time=30:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=../98_logs/ana_%j.out
#SBATCH --error=../98_logs/ana_%j.err

../../software/bedtools2/bin/bedtools \
	intersect \
	-a ../18_ppi_interface/ukb_rand.bed \
	-b ../18_ppi_interface/ppi_interface.bed \
	-c > ../18_ppi_interface/ukb_rand_at_interface.bed


