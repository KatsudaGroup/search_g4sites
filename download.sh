#/bin/sh

organ_list=(human mouse)
data_dir=data

if [ ! -e $data_dir ]; then
	mkdir $data_dir
fi

if printf '%s\n' "${organ_list[@]}" | grep -qx "human"; then
	for i in `seq 1 12`;
	do
		filename=human.${i}.rna.fna
		url=https://ddbj.nig.ac.jp/public/mirror_database/refseq/H_sapiens/mRNA_Prot/${filename}.gz
		
		if [ ! -e ${data_dir}/${filename} ]; then
			wget -P ${data_dir} ${url}
			gunzip ${data_dir}/${filename}.gz
		fi
	done;
fi
