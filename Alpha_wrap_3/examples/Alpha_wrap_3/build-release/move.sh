model=$1
omega=$2

results=results
mkdir -p $results

expe=$results/$model-$omega
rm -rf $expe
mkdir $expe

data=$expe/data
mkdir $data
mv *xyz *poly *off $data
