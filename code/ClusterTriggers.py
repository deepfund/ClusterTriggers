import sys, os.path, argparse
from math import log, exp
import numpy as np
import configparser

def dispersion_sec_ghz(v0, v1, dm):
    return abs(4.15 * 1E-3 * dm * (v0 ** -2 - v1 ** -2))


def make_default_config():

    config = configparser.ConfigParser()

    config['io'] = {
        'input_file': '',
        'output_file': ''
    }
    
    config['dm'] = {
        'dm_axis_min': '0.1',
        'dm_axis_max': '1000',
        'dm_axis_count': '128',
        'dm_axis_log': 'yes'
    }
    
    config['time'] = {
        'time_step': '0.1',
        'max_cluster_duration': '0.1',
        'time_buffer_size': '0.5'
    }

    
    with open('default.ini', 'w') as configfile:
        config.write(configfile)

# -------------------------------------------------------------------------------------
# A helper class to convert DM values to discrete grid locations
# usage:
#   dmc = DiscretizeDm(args)
#   dmi = dmc.convert(40.2)
# -------------------------------------------------------------------------------------
class DiscretizeDm:

    def __init__(self, args):
        self.dm_min = args.dm_min
        self.dm_max = args.dm_max
        self.dm_steps = args.dm_steps
        self.dm_log = args.dm_log
        
        if self.dm_log:
            self.convert = self.convert_log
            self.dm_min = max(0.01, self.dm_min)        # a hard lower bound
            self.log_dm_min = log(self.dm_min)
            self.log_dm_max = log(self.dm_max)
            self.dm_delta = (self.log_dm_max - self.log_dm_min) / self.dm_steps
        else:
            self.convert = self.convert_linear
            self.dm_delta = (self.dm_max - self.dm_min) / self.dm_steps

    def convert_linear(self, dm):
        dmi = (dm - self.dm_min) / self.dm_delta
        dmi = int(dmi)
        dmi = max(0, dmi)
        dmi = min(dmi, self.dm_steps - 1)
        return dmi
        
    def convert_log(self, dm):
        dmi = (log(dm) - self.log_dm_min) / self.dm_delta
        dmi = int(dmi)
        dmi = max(0, dmi)
        dmi = min(dmi, self.dm_steps - 1)
        return dmi


# -------------------------------------------------------------------------------------
# Cluster engine
# -------------------------------------------------------------------------------------
class ClusterEngine:

    def __init__(self, integration_step, args):
        self.debug = False
        self.integration_step = integration_step
        self.args = args
        self.time_step = args.time_step
        self.dm_steps = args.dm_steps
        self.buffer_size = args.buffer_size

        self.discretize_dm = DiscretizeDm(args)
        self.clusters = {}
        self.max_cluster_id = 0

        # setup the grid
        # We discretize {t,dm} values onto a 2D grid for fast lookup and clustering.
        # the grid contains cluster id, each grid cell can be empty (value = 0) or covered by a cluster with that id
        # The grid is a sliding window, we write solely in the right (R) section of the grid
        # if new data falls outside the grid on the right then we shift the grid with a half width nr units until the
        # the data would fall again on the grid.
        # If new data falls outside the grid on the left we drop the trigger..
        # We add an extra empty row and column to each side of the grid to simplify the neighborhood search code.
        #
        # Below an example of a grid row with a 2x 4 elements (the buffer size)
        #   0   1 2 3 4   5 6 7 8   9
        # [ x | L L L L | R R R R | x ]
        # grid_half_width = 4, grid_width = 10

        self.grid_height = self.dm_steps + 2

        # The half width is 4 in the example above
        self.grid_half_width = 1 + int(self.buffer_size / self.time_step)

        # The full width including an extra row om both sides is 10 in the example above
        self.grid_width = 2 * self.grid_half_width + 2

        # set the initial write head to pos 5 in the example above.
        self.grid_t0 = 0

        # Allocate memory for the grid
        self.grid = np.zeros((self.grid_width, self.grid_height), dtype=np.int)
        if self.debug:
            print('NEW ClusterEngine {} {}'.format(
                self.integration_step,
                self.args
            ))

    def shift_grid(self, shift_steps, pos_i):
        if self.debug:
            print('shifting grid', shift_steps, 'because pos_i', pos_i, 'grid_width', self.grid_width)

        if shift_steps == 0:
            return pos_i

        if shift_steps == 1:
            # copy the right half to the left half
            self.grid[1:self.grid_half_width + 1] = self.grid[self.grid_half_width + 1:-1]

            # clear the old right half
            self.grid[self.grid_half_width + 1:-1] = 0

        elif shift_steps > 1:
            self.grid[:, :] = 0

        elif shift_steps == -1:
            self.grid[self.grid_half_width + 1:-1] = self.grid[1:self.grid_half_width + 1]

        elif shift_steps < -1:
            self.grid[:, :] = 0

        # move the origin
        self.grid_t0 += shift_steps * self.grid_half_width
        pos_i -= shift_steps * self.grid_half_width
        if self.debug:
            print('after the shift pos_i is', pos_i)
        return pos_i

    # if our x (time) position on the grid is outside the bound then we need to shift the grid
    def calc_grid_shift(self, pos_i):
        if pos_i >= self.grid_width - 1:  # are we past the right on the grid?
            return 1 + int((1 + pos_i - self.grid_width) / self.grid_half_width)
        if pos_i <= 0:  # are we past the left on the grid?
            return -1 - int(-pos_i / self.grid_half_width)
        return 0

    # This function discetizes time to positions on the grid, and optionally move the grid to make sure
    # that this even will fall on the right part of the grid.
    # Return the position of t on the the grid.
    def compute_grid_x(self, t):
        # Discretize time
        t_i = int(t / self.time_step)
        pos_i = t_i - self.grid_t0

        # Check if we need to shift the grid
        shift_steps = self.calc_grid_shift(pos_i)
        if shift_steps != 0:
            pos_i = self.shift_grid(shift_steps, pos_i)

        return pos_i

    def handle(self, t, dm, snr, sample):
        
        # compute the position of this event on the discrete grid
        grid_x = self.compute_grid_x(t)
        #grid_y = int(min(dm, self.dm_max) / self.dm_step) + 1
        grid_y = self.discretize_dm.convert(dm)

        if self.debug:
            print('handle', grid_x, grid_y, t, dm, sample)

        # check if this grid cell is part of a cluster
        c1 = self.grid[grid_x, grid_y]
        if c1 > 0:
            self.join_cluster(c1, t, dm, snr, sample)
            return

        # find the ids of the neighboring clusters
        nn_subgrid = self.grid[grid_x - 1:grid_x + 2, grid_y - 1:grid_y + 2]
        nn_cluster_ids = np.unique(nn_subgrid[nn_subgrid > 0])

        #  there is no cluster close-by: we create a new one
        if len(nn_cluster_ids) == 0:
            self.new_cluster(t, dm, snr, sample)
            self.grid[grid_x, grid_y] = self.max_cluster_id
            return

        c1 = nn_cluster_ids[0]
        # merge cluster if there are multiple neightbor clusters
        if len(nn_cluster_ids) > 1:
            for j in range(1, len(nn_cluster_ids)):
                c2 = nn_cluster_ids[j]
                self.merge_custers(c1, c2)

        # now join the main cluster
        self.join_cluster(c1, t, dm, snr, sample)
        self.grid[grid_x, grid_y] = c1

    @staticmethod
    def print_cluster_header():
        print('integration_step t_min t_max dm_min dm_max count peak_snr peak_t peak_dm peak_sample')

    def print_cluster(self, c1):
        print('{} {:.4f} {:.4f} {:.4f} {:.4f} {} {:.4f} {:.4f} {:.4f} {}'.format(
            self.integration_step,
            self.clusters[c1]['t_min'],
            self.clusters[c1]['t_max'],
            self.clusters[c1]['dm_min'],
            self.clusters[c1]['dm_max'],
            self.clusters[c1]['count'],
            self.clusters[c1]['peak_snr'],
            self.clusters[c1]['peak_t'],
            self.clusters[c1]['peak_dm'],
            self.clusters[c1]['peak_sample']
        ))

    def new_cluster(self, t, dm, snr, sample):
        self.max_cluster_id += 1
        self.clusters[self.max_cluster_id] = {
            'id': self.max_cluster_id, 
            't_min': t, 
            't_max': t, 
            'dm_min': dm,
            'dm_max': dm, 
            'peak_snr': snr, 
            'peak_t': t, 
            'peak_dm': dm, 
            'count': 1, 
            'peak_sample': sample
        }
        if self.debug:
            print('created a new cluster with id', self.max_cluster_id, 'and joined it')

    def join_cluster(self, c1, t, dm, snr, sample):
        self.clusters[c1]['t_min'] = min(self.clusters[c1]['t_min'], t)
        self.clusters[c1]['t_max'] = max(self.clusters[c1]['t_max'], t)
        self.clusters[c1]['dm_min'] = min(self.clusters[c1]['dm_min'], dm)
        self.clusters[c1]['dm_max'] = max(self.clusters[c1]['dm_max'], dm)
        self.clusters[c1]['count'] = self.clusters[c1]['count'] + 1
        if snr > self.clusters[c1]['peak_snr']:
            self.clusters[c1]['peak_snr'] = snr
            self.clusters[c1]['peak_t'] = t
            self.clusters[c1]['peak_dm'] = dm
            self.clusters[c1]['peak_sample'] = sample
        if self.debug:
            print('joined cluster with id', c1)

    def merge_custers(self, c1, c2):

        self.clusters[c1]['t_min'] = min(self.clusters[c1]['t_min'], self.clusters[c2]['t_min'])
        self.clusters[c1]['t_max'] = max(self.clusters[c1]['t_max'], self.clusters[c2]['t_max'])
        self.clusters[c1]['dm_min'] = min(self.clusters[c1]['dm_min'], self.clusters[c2]['dm_min'])
        self.clusters[c1]['dm_max'] = max(self.clusters[c1]['dm_max'], self.clusters[c2]['dm_max'])
        
        if self.clusters[c1]['peak_snr'] < self.clusters[c2]['peak_snr']:
            self.clusters[c1]['peak_snr'] = self.clusters[c2]['peak_snr']
            self.clusters[c1]['peak_t'] = self.clusters[c2]['peak_t']
            self.clusters[c1]['peak_dm'] = self.clusters[c2]['peak_dm']
            self.clusters[c1]['peak_sample'] = self.clusters[c2]['peak_sample']

        # update the cluster ids on the grid
        self.grid[self.grid == c2] = c1

        # remove the old cluster info
        self.clusters.pop(c2, None)
        if self.debug:
            print('merge_clusters', c1, c2)

# -------------------------------------------------------------------------------------
# The 5 col format parser:
# "# DM      Sigma      Time (s)     Sample    Downfact"
# -------------------------------------------------------------------------------------
def parse_format_C5(line):
    va = line.strip().split()
    beam = 0
    batch = 0
    sample = va[3]
    integration_step = va[4]
    t = va[2]
    DM = va[0]
    SNR = va[1]
    return beam, batch, sample, integration_step, t, DM, SNR

# -------------------------------------------------------------------------------------
# The 7 col format parser:
# " # beam batch sample integration_step time DM SNR
# -------------------------------------------------------------------------------------
def parse_format_C7(line):
    va = line.strip().split()
    beam = int(va[0])
    batch = int(va[1])
    sample = int(va[2])
    integration_step = int(float(va[3]))
    t = float(va[4])
    DM = float(va[5])
    SNR = float(va[6])
    return beam, batch, sample, integration_step, t, DM, SNR

# -------------------------------------------------------------------------------------
# The 10 col format parser:
# "# beam_id batch_id sample_id integration_step compacted_integration_steps time DM_id DM compacted_DMs SNR"
# -------------------------------------------------------------------------------------
def parse_format_C10(line):
    va = line.strip().split()
    beam = int(va[0])
    batch = int(va[1])
    sample = int(va[2])
    integration_step = int(float(va[3]))
    t = float(va[5])
    DM = float(va[7])
    SNR = float(va[9])
    return beam, batch, sample, integration_step, t, DM, SNR


# -------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# -------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(
        description='This tool clusters AMBER triggers',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-i', '--infile',
        nargs='?',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='Filename to read triggers from, reads from stdin when omitted.'
    )

    parser.add_argument(
        '-o', '--outfile',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Filename to write clusters to, writes to stdout when omitted.'
    )

    parser.add_argument(
        '-ts', '--time-step',
        default=0.1,
        type=float,
        help=('Time discretisation step (seconds). During clustering all events timestamps are rounded '
              'to multiples of this value.')
    )

    parser.add_argument(
        '-dms', '--dm-steps',
        default=64,
        type=int,
        help=('Number of dispersion measure discretisation step. During clustering all events '
        'dispersion measures are rounded to multiples of this value.')
    )
    
    parser.add_argument(
        '-dmlo', '--dm-min',
        default=1,
        type=float,
        help='Maximum dispersion measure, all values above will be treated as if there where capped at this value'
    )

    parser.add_argument(
        '-dmhi', '--dm-max',
        default=1000,
        type=float,
        help='Maximum dispersion measure, all values above will be treated as if there where capped at this value'
    )

    parser.add_argument(
        '-dmlog', '--dm-log',
        type=str2bool,
        nargs='?',
        const=True,
        default=True,
        help="Use log scale for the DM axis."
    )
    parser.add_argument(
        '-ct', '--max-duration',
        default=0.2,
        type=float,
        help=('Maximum time duration of a cluster (seconds). If a cluster tries to grow bigger than this timespan '
              'the cluster will be broken in multiple clusters.')
    )

    #parser.add_argument(
    #    '-z', '--buffer-size',
    #    default=5,
    #    type=float,
    #    help=('Maximum time (seconds) that we hang on to cluster  for modification before we flush then to '
    #          'the output file.')
    #)

    parser.add_argument(
        '-b', '--beam',
        type=int,
        help='Only process triggers from a certain beam.'
    )

    parser.add_argument(
        '-is', '--integration-step',
        type=int,
        help='Only process triggers that have a certain integration step.'
    )

    parser.add_argument(
        '-iis', '--ignore-integration-step',
        type=str2bool,
        nargs='?',
        const=True,
        default=True,
        help='Ignore integration step information.'
    )

    parser.add_argument(
        '-siglo', '--sig-thresh',
        default=5.0,
        type=float,
        help='Minimal SNR'
    )
    parser.add_argument(
        '-sighi', '--sig-max',
        default=np.inf,
        type=float,
        help='Maximum SNR'
    )

    parser.add_argument(
        '-mr', '--max-rows',
        default=np.inf,
        type=int,
        help='Maximum number of rows to read'
    )
    
    args = parser.parse_args()
    args.buffer_size = 5
    return args
# -------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------
def main(args):

    # dictionary of cluster engine, one for each integration_step
    cluster_engines = {}

    # ---------------------------------------------------------------------------------
    # Determine input file format
    # ---------------------------------------------------------------------------------
    parser = None
    hdr = args.infile.readline().strip().split()
    
    if not isinstance(hdr, (list,)):
        print("Error determining trigger file format.")
        print("Expected space delimited column names.")
        exit(1)

    if len(hdr)==0 or hdr[0] != '#':
        print("Error determining trigger file format.")
        print("The first row of he trigger file doesn't start with a #. ")
        exit(1)
    
    if hdr == ['#', 'DM', 'Sigma', 'Time', '(s)', 'Sample', 'Downfact']:
        parser = parse_format_C5
    
    if hdr == ['#', 'beam', 'batch', 'sample', 'integration_step', 'time', 'DM', 'SNR']:
        parser = parse_format_C7
    
    if hdr == ['#', 'beam_id', 'batch_id', 'sample_id', 'integration_step', 'compacted_integration_steps', 'time', 'DM_id', 'DM', 'compacted_DMs', 'SNR']:
        parser = parse_format_C10
    
    if not parser:
        print("Error determining trigger file format.")
        print("Unknown header format.")
        print(hdr)
        exit(1)
    
    # process trigger lines
    row_count = 0
    for line in args.infile:

        # check if we have a row count limit
        row_count += 1
        if row_count > args.max_rows:
            break
            
        try:
            beam, batch, sample, integration_step, t, DM, SNR = parser(line)
            
            

            # -------------------------------------------------------------
            # apply filters, discard trigger
            # -------------------------------------------------------------
            if args.integration_step is not None:
                if int(args.integration_step) != integration_step:
                    continue
            
            if args.beam is not None:
                if int(args.beam) != beam:
                    continue
            
            if args.sig_thresh is not None:
                if SNR < args.sig_thresh:
                    continue
            
            if args.sig_max is not None:
                if SNR > args.sig_max:
                    continue

            # -------------------------------------------------------------
            # use 1 or multiple engines for clustering
            # -------------------------------------------------------------
            if args.ignore_integration_step:
                engine_name = 'generic'
            else:
                engine_name = 'is_{}'.format(integration_step)

            # -------------------------------------------------------------
            # check if we already have a cluster engine for this integration step
            # -------------------------------------------------------------
            if not engine_name in cluster_engines:
                ceng = ClusterEngine(integration_step, args)
                cluster_engines[engine_name] = ceng
            cluster_engines[engine_name].handle(t, DM, SNR, sample)
        except:
            pass

    
    # Print the clusters
    if not args.outfile or args.outfile == sys.stdout:
        ClusterEngine.print_cluster_header()
        for k, v in cluster_engines.items():
            for k2, v2 in v.clusters.items():
                v.print_cluster(v2['id'])

    sig_cut, dm_cut, tt_cut, ds_cut, ind_full = [],[],[],[],[]
    for k, v in cluster_engines.items():
        for k2, v2 in v.clusters.items():
        
            # ---------------------------------
            # available cluster info:
            # ---------------------------------
            # v.integration_step
            # v2['id']
            # v2['t_minv
            # v2['t_max']
            # v2['dm_min']
            # v2['dm_maxv
            # v2['count']
            # v2['peak_snr']
            # v2['peak_t']
            # v2['peak_dm']
            # v2['peak_sample']
            
            sig_cut.append(v2['peak_snr'])
            dm_cut.append(v2['peak_dm'])
            tt_cut.append(v2['peak_t'])
            ds_cut.append(v.integration_step)
            ind_full.append(v2['peak_sample'])
    
    
    return np.array(sig_cut), np.array(dm_cut), np.array(tt_cut), np.array(ds_cut), np.array(ind_full)

# ----------------------------------------------------------------------------------------
# Interface
# ----------------------------------------------------------------------------------------
def get_triggers(
        infile, 
        sig_thresh=5.0, 
        dm_min=0, 
        dm_max=np.inf, 
        t_window=0.5, 
        max_rows=None, 
        t_max=np.inf,
        sig_max=np.inf,
        dt = 40.96,
        delta_nu_MHz=300./1536, 
        nu_GHz=1.4,
        fnout=False
    ):
    
    args = get_args()
    
    args.infile = open(infile, 'r')
    args.dm_min = dm_min
    args.dm_max = dm_max
    args.sig_max = sig_max
    args.sig_thresh = sig_thresh
    args.max_rows = np.inf
    
    sig_cut, dm_cut, tt_cut, ds_cut, ind_full = main(args)
    
    if fnout != False:
        clustered_arr = np.concatenate([sig_cut, dm_cut, tt_cut, ds_cut, ind_full])
        clustered_arr = clustered_arr.reshape(5, -1)
        np.savetxt(fnout, clustered_arr) 
        
    return sig_cut, dm_cut, tt_cut, ds_cut, ind_full

# ----------------------------------------------------------------------------------------
# Command line version
# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = get_args()
    main(args)
