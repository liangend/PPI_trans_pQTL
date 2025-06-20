source ~/.bashrc
conda activate ldsc

sumstats=$1
ref_annot=$2
weight_ld=$3
frqfile=$4
h2_output=$5

python /project/xuanyao/jinghui/software/ldsc/ldsc.py \
	--h2 $sumstats \
	--ref-ld-chr $ref_annot \
	--w-ld-chr $weight_ld \
	--overlap-annot \
	--frqfile-chr $frqfile \
	--print-coefficients \
	--no-intercept \
	--out $h2_output

