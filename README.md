# ClusterTriggers
Amber cluster triggering tools


Usage:

for reading/writing to files:

````
python3 ClusterTriggers.py --infile test.trigger --integration-step 25 --outfile clusters.txt
````

for I/O streams
````
cat test.trigger | python3 ClusterTriggers.py >> clusters.txt
````


IN the python shell:

````
import ClusterTriggers

fntrig = 'sample.trigger'
snr, dm, time_arrival, downsample, indexes = ClusterTriggers.get_triggers(fntrig, sig_thresh=8.0, max_rows=1000000, fnout='sample.clusters')
````
