#!/bin/bash

# use traditional way for conda environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate UniFold


#setting environment for cuda-11.0 gpu2-5
export PATH=/usr/local/cuda-11.4/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH

# User configuration
db_dir=/mnt/db;
pretrained_data_dir=/mnt/db/unifold/

# automatically determined directory
af_official_repo=/repo/alphafold/ ;
uf_official_repo=$(readlink -f $(dirname $0)) ;
dir=`pwd`;


out_dir=$dir/output;
res_dir=$dir/res;



usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        # edited by Yinying
        echo "-m <model_preset>  Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model"
        echo "-n <num_multimer_predictions_per_model>       How many predictions (each with a different random seed) will be generated per model"
        echo "-t <template_date> Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo "-p <pretrained_data_date> Pretrained data release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo "-r <run_relax> Pretrained data release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo ""
        exit 1
}

while getopts ":m:t:n:e:p:r:" i; do
        case "${i}" in

        m)
                model_preset=$OPTARG
        ;;

        t)
                max_template_date=$OPTARG
        ;;
        n)
                num_multimer_predictions_per_model=$OPTARG
        ;;
        e)
                num_ensemble=$OPTARG
        ;;
        p)
                pretrained_data_date=$OPTARG
        ;;
        r)
                run_relax=$OPTARG
        ;;
        *)
                echo Unknown argument!
                usage
        ;;

        esac
done

if [[ "$max_template_date" == "" ]] ; then
    max_template_date=2021-10-30
fi

if [[ "$max_template_date" == "no" ]] ; then
    max_template_date=1900-01-01
fi


if [[ "$pretrained_data_date" == ""  ]] ; then
    pretrained_data_date=2022-08-01
elif [[ ! -d ${pretrained_data_dir}/${pretrained_data_date}/ ]];then
            echo "ERROR: Unknown pretrained_data_date ${pretrained_data_date} or the pretrained_data_dir ${pretrained_data_dir} inaccessible. "
            usage

fi



# edited by Yinying
if [[ "$model_preset" == "" ]] ; then
    model_preset="monomer"
fi

# set default num_ensemble
if [[ "$model_preset" == "monomer" || "$model_preset" == "monomer_ptm" || "$model_preset" == "multimer" ]] ; then
    if [[ "$num_ensemble" == "" ]] ; then
        num_ensemble=1
    fi
fi

if [[ "$model_preset" == "monomer_casp14"  ]] ; then
    if [[ "$num_ensemble" == "" ]] ; then
        num_ensemble=8
    fi
fi

if [[ "$model_preset" == "multimer" ]] ; then
    if [[ "$num_multimer_predictions_per_model" == "" ]];then
        num_multimer_predictions_per_model=2
    fi
else
    num_multimer_predictions_per_model=1
fi

if [[ "$run_relax" == "" || "$run_relax" == "true" ]] ; then
    run_relax=true
else
    run_relax=false
fi

if [[ "$model_preset" != "monomer" && "$model_preset" != "monomer_casp14" && "$model_preset" != "monomer_ptm" && "$model_preset" != "multimer" ]] ; then
    echo "Unknown model_preset! "
    usage
fi

# preset configurations
#fasta_path=$1
#output_dir_base=$2
database_dir=$db_dir
max_template_date=$max_template_date
param_path=$pretrained_data_dir



echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"



mkdir $res_dir;
mkdir $out_dir;
mkdir $res_dir/models;
mkdir $res_dir/lite;
mkdir $res_dir/full;
mkdir $dir/processed;


UF_process(){
	local dir=$1;
	local i=$2;
	local decoy_name=${i%.fasta};

	#cd $af_official_repo; # fix error caused by some path configs.
	pushd $out_dir
    if [ ! -f "$res_dir/lite/${decoy_name}_UF2_lite.tar.bz2" ]; then
		# not started yet??
		if [ $(ls  $out_dir/$decoy_name/UF | grep -e ".feature.pkl.gz$" |wc -l) -lt 1 ]; then
            echo "File does not exist: $out_dir/$decoy_name/UF/*.feature.pkl.gz";
        	echo Modeling is not started: $i;
        	echo "Starting homogeneous searching..."
	        local cmd="python ${uf_official_repo}/unifold/homo_search.py \
              --fasta_path=$dir/$i \
              --max_template_date=$max_template_date \
              --output_dir=$decoy_name  \
              --uniref90_database_path=$database_dir/uniref90/uniref90.fasta \
              --mgnify_database_path=$database_dir/mgnify/mgy_clusters.fa \
              --bfd_database_path=$database_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
              --uniclust30_database_path=$database_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
              --uniprot_database_path=$database_dir/uniprot/uniprot.fasta \
              --pdb_seqres_database_path=$database_dir/pdb_seqres/pdb_seqres.txt \
              --template_mmcif_dir=$database_dir/pdb_mmcif/mmcif_files \
              --obsolete_pdbs_path=$database_dir/pdb_mmcif/obsolete.dat \
              --use_precomputed_msas=True";
	        echo "$cmd";eval "$cmd"
        else
            echo Find feature files in $out_dir/$decoy_name/features.pkl;
            echo Skip the run_feature process: $decoy_name
        fi
        # featuring is not started yet, Skip this sequence.
		if [ $(ls  $out_dir/$decoy_name/UF | grep -e ".feature.pkl.gz$" |wc -l) -lt 1 ]; then
            echo "Feature files not found: $out_dir/$decoy_name/UF/*.feature.pkl.gz"
        	echo Modeling is not started: $i;
        	echo Skip Modeling: $i;
	    else
	        #
            echo "Find feature files: $out_dir/$decoy_name/UF/*.feature.pkl.gz";
            echo Runing modeling process: $decoy_name
            if [[ -f "${pretrained_data_dir}/${model_preset}.unifold.pt" ]];then
              echo "Pretrained model file not found: ${pretrained_data_dir}/${model_preset}.unifold.pt";
              exit 1;
            fi

            if [ $(ls  $out_dir/$decoy_name/UF | grep -e ".pdb$" |wc -l) -lt 1 ]; then
              cmd="python ${uf_official_repo}/unifold/inference.py \
                --param_path=${pretrained_data_dir}/${pretrained_data_date}/${model_preset}.unifold.pt \
                --data_dir=${decoy_name} \
                --model_name=multimer_ft \
                --target_name=UF \
                --output_dir=${decoy_name}" ;

              echo "$cmd";eval "$cmd"
            fi

		    cd $out_dir && \
		    echo Selecting best_model ...
		    pushd $decoy_name/UF
		        for f in *.pdb;do cp ${f} ${res_dir}/models/${decoy_name}_${f};done
		    popd

		    # plot summary of af modeling results
		    mkdir -p $res_dir/models/plot/$decoy_name
		    python ${af_official_repo}/AlphaPickle.py --res_dir $decoy_name/UF --save_dir $res_dir/models/plot/$decoy_name

		    #cp $decoy_name/ranked_0.pdb $res_dir/best_model/${decoy_name}_ranked_0.pdb &&
#		    echo Collecting results files .... && \
#		    tar jcf $decoy_name\_UF2_lite.tar.bz2  --exclude *.pkl.gz --exclude $decoy_name/msas $decoy_name && \
#		    mv $decoy_name\_UF2_lite.tar.bz2 $res_dir/lite && \
#		    tar jcf $decoy_name\_UF2_full.tar.bz2 --remove-files $decoy_name &&  \
#		    mv $decoy_name\_UF2_full.tar.bz2 $res_dir/full && \
#		    mv $dir/$i $dir/processed &
        fi
    else
            echo "Find modeling files in $res_dir/lite/${decoy_name}_UF2_lite.tar.bz2";
            echo Skip the run_alphafold process: $decoy_name
    fi
    popd
}

cd $dir;
echo We are now in `pwd`

total=`ls  |grep .fasta |wc -l`;
fin=0;
rest=$total;

for i in `ls  |grep .fasta`;
do
    echo $dir/$i;
    cd $dir;
    # run processing ...
    UF_process $dir $i ;
    # mv $dir/$i $dir/processed;
    let fin++;
    let rest--;

done

wait;

echo Sending final notify ....
echo python ${af_official_repo}/sms.py `whoami` `basename $dir`  $total $fin $rest;

echo "++++++++++++++++++++++++++++++++++++++++"

echo "process end at : "
date
cd $dir;

