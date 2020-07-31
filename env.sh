module load boost/1.66.0
module load python/3.7

function theia(){
  source /global/scratch/users/maskins/snoplus/env.sh
  source /global/scratch/users/maskins/WblsPathBTools/ratpac/env.sh
  export PATH=/global/scratch/users/maskins/WblsPathBTools/off_recon_nlopt:/global/scratch/users/maskins/WblsPathBTools/bin:$PATH
  source /global/scratch/users/maskins/WblsPathBTools/off_recon_nlopt/env-recon.sh
}

theia
