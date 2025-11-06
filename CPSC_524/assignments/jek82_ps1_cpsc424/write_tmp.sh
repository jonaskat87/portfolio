# run bash script once Slurm job has been completed; it will write output of triad into tmp.txt for plotting in MATLAB
cd /home/cpsc424_jek82/project/assignments/jek82_ps1_cpsc424
# get line number where vector triad benchmarking starts
Ni=`awk '{if ($1=="***Results" && $3=="vector") print NR}' slurm-ps1.out-*`
# get line number where CPU architecture starts
Nf=`awk '{if ($1=="***CPU") print NR}' slurm-ps1.out-*`
awk -v start="$Ni" -v end="$Nf" '{if (NR>start && NR<end-1) print}' slurm-ps1.out-* > tmp.txt
