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

# Command Line arguments
```
usage: ClusterTriggers.py [-h] [-i [INFILE]] [-o [OUTFILE]] [-ts TIME_STEP]
                          [-dms DM_STEPS] [-dmlo DM_MIN] [-dmhi DM_MAX]
                          [-dmlog [DM_LOG]] [-ct MAX_DURATION]
                          [-z BUFFER_SIZE] [-b BEAM] [-is INTEGRATION_STEP]
                          [-iis [IGNORE_INTEGRATION_STEP]] [-siglo SIG_THRESH]
                          [-sighi SIG_MAX] [-mr MAX_ROWS]


optional arguments:
  -h, --help            show this help message and exit
  -i [INFILE], --infile [INFILE]
                        Filename to read triggers from, reads from stdin when
                        omitted. (default: <_io.TextIOWrapper name='<stdin>'
                        mode='r' encoding='UTF-8'>)
  -o [OUTFILE], --outfile [OUTFILE]
                        Filename to write clusters to, writes to stdout when
                        omitted. (default: <_io.TextIOWrapper name='<stdout>'
                        mode='w' encoding='UTF-8'>)
  -ts TIME_STEP, --time-step TIME_STEP
                        Time discretisation step (seconds). During clustering
                        all events timestamps are rounded to multiples of this
                        value. (default: 0.1)
  -dms DM_STEPS, --dm-steps DM_STEPS
                        Number of dispersion measure discretisation step.
                        During clustering all events dispersion measures are
                        rounded to multiples of this value. (default: 64)
  -dmlo DM_MIN, --dm-min DM_MIN
                        Maximum dispersion measure, all values above will be
                        treated as if there where capped at this value
                        (default: 1)
  -dmhi DM_MAX, --dm-max DM_MAX
                        Maximum dispersion measure, all values above will be
                        treated as if there where capped at this value
                        (default: 1000)
  -dmlog [DM_LOG], --dm-log [DM_LOG]
                        Use log scale for the DM axis. (default: True)
  -ct MAX_DURATION, --max-duration MAX_DURATION
                        Maximum time duration of a cluster (seconds). If a
                        cluster tries to grow bigger than this timespan the
                        cluster will be broken in multiple clusters. (default:
                        0.2)
  -b BEAM, --beam BEAM  Only process triggers from a certain beam. (default:
                        No filtering, all beams)
  -is INTEGRATION_STEP, --integration-step INTEGRATION_STEP
                        Only process triggers that have a certain integration
                        step. (default: No filtering, all integration step)
  -iis [IGNORE_INTEGRATION_STEP], --ignore-integration-step [IGNORE_INTEGRATION_STEP]
                        Ignore integration step information. (default: True)
  -siglo SIG_THRESH, --sig-thresh SIG_THRESH
                        Minimal SNR (default: 5.0)
  -sighi SIG_MAX, --sig-max SIG_MAX
                        Maximum SNR (default: inf)
  -mr MAX_ROWS, --max-rows MAX_ROWS
                        Maximum number of rows to read (default: inf)
