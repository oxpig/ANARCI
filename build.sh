clean ()
{
	rm -rf anarci.egg-info/
	rm -rf dist/
	rm -rf build/
	rm -rf build_pipeline/IMGT_sequence_files/
	rm -rf build_pipeline/muscle_alignments/
	rm -rf build_pipeline/curated_alignments/
	rm -rf build_pipeline/HMMs/

	rm -f lib/python/anarci/germlines.py
	rm -rf lib/python/anarci/dat
}

if ! [ -z $1 ]; then
	if [ $1 = "clean" ]; then
		clean
	else
		echo "Unknown action. Use no argument to build and 'clean' to remove built files. Exiting."
	fi

	exit 0
fi

cd build_pipeline/
echo "Downloading germlines from IMGT and building HMMs..."
bash RUN_pipeline.sh

cp curated_alignments/germlines.py ../lib/python/anarci/
mkdir -p ../lib/python/anarci/dat
cp -r HMMs ../lib/python/anarci/dat

