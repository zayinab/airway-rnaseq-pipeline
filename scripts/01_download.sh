# Download the first 1 million read pairs of each airway sample
cd data/raw
for id in $(cat ../samples.txt); do
  echo "Downloading $id ..."
  fastq-dump -X 1000000 --split-files --gzip $id
done

