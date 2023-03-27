# Read Input
# Get the batch number
batch_num=$1
shift

# Get working directory
CWD=$(pwd)
echo "Working directory: $CWD"

# Make Directories
echo "Creating Directories For: PBS Scripts/Outfiles, Intermediate ChooseR data files, ChooseR results"
mkdir -p batch_pbs_$batch_num
mkdir -p batch_data_$batch_num
mkdir -p batch_results_$batch_num

# Generate a pbs script
generate_pbs() {
	min_dist=$1
	n_neighbors=$2
	min_dist_str=${min_dist//./}
	#n_neighbors=${n_neighbors//./_}

	pbs_str="pbs_$min_dist_str""_$n_neighbors"".sbatch"
	echo $pbs_str
	echo "#!/bin/bash
#SBATCH -A gts-js585-biocluster
#SBATCH -Jclust_$min_dist_str""_$n_neighbors
#SBATCH --mem=128gb
#SBATCH -N 2 --ntasks-per-node=4
#SBATCH -t 44:00:00
#SBATCH -qinferno
#SBATCH -opbs_$min_dist_str""_$n_neighbors"".out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ggruenhagen3@gatech.edu

cd \$SLURM_SUBMIT_DIR
source ~/.bashrc

local_dir=~/scratch/brain/chooser/chooser_spatial/examples/

Rscript \$local_dir/hb_chooser.R $batch_num $min_dist $n_neighbors 5" > $pbs_str
sbatch $pbs_str
}

# Clustering Parameters to Loop Through
# min_dist parameters
declare -a min_dist_array=("0.001" "0.1" "0.2" "0.3" "0.4" "0.5" )
# n_neighbors parameters
declare -a n_neighbors_array=("5" "10" "20" "30" "40" "50" )

# Generate pbs scripts
cd batch_pbs_$batch_num 
for this_min_dist in ${min_dist_array[@]}; do
	for this_n_neighbors in ${n_neighbors_array[@]}; do
		generate_pbs $this_min_dist $this_n_neighbors
	done
done
cd $CWD

