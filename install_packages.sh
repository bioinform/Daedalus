##first download packages.tar.gz
##then 'tar -xvfz packages.tar.gz'
cd ./packages/fastq-streamer && python setup.py install && cd ../../
cd ./packages/bam-streamer && python setup.py install && cd ../../
cd ./packages/trim-primers && python setup.py install && cd ../../ 
cd ./packages/parse-umi && python setup.py install && cd ../../ 
cd ./packages/SeqNetworks && python setup.py install && cd ../../ 
cd ./packages/ipete-dedup && python setup.py install && cd ../../
cd ./packages/extract-umi && python setup.py install && cd ../../
cd ./packages/ipete-metrics && python setup.py install && cd ../../ 
cd ./packages/ipete-reporter && python setup.py install && cd ../../ 
cd ./packages/ipete-spikein-parser && python setup.py install && cd ../../ 
cd ./packages/spikein-split && python setup.py install && cd ../../
cd ./packages/VDJ-recombinant-detection && python setup.py install && cd ../../

