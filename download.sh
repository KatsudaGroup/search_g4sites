#/bin/sh
#
data_dir=data
if [ ! -e $data_dir ]; then
	mkdir $data_dir
fi

#wget -P ${data_dir} https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.files.installed
wget -P ${data_dir} https://ddbj.nig.ac.jp/public/mirror_database/refseq/H_sapiens/mRNA_Prot/human.files.installed

filename_pattern="human.*.rna.fna.gz"
IFS=$'\t'
while read -r md5_value filename; do
	#echo ${filename}
	case "${filename}" in 
		${filename_pattern})
			echo ${filename}
			#url=https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/${filename}
			url=https://ddbj.nig.ac.jp/public/mirror_database/refseq/H_sapiens/mRNA_Prot/${filename}
			echo "download ${filename}"
			wget -P ${data_dir} ${url}

			actual_hash=`md5 -q ${data_dir}/${filename}`
			echo $actual_hash
			if [[ $actual_hash == $md5_value ]]; then
				echo "MD5 Match: ${filename}"
				gunzip ${data_dir}/${filename}
			else
				echo "Download but MD5 Did not matched!"
				break
			fi
			;;

		*)
			;;
	esac

done < ${data_dir}/human.files.installed

