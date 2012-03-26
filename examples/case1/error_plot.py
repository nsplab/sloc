#!/usr/bin/env python
"""
Plot source localization errors vs SNR
"""

import case1
import numpy
import pylab

dist = lambda x,y: sum((x-y)**2)**0.5

def read_dipoles(filename):
    dipoles = []
    with open(filename, 'r') as fp:
        n = int(fp.readline())
        for i in range(n):
            dipoles.append(numpy.array([float(x) for x in fp.readline().split()]))
    return dipoles


class ErrorTable(object):

    def __init__(self):
        self.dipoles = dict()
        self.sources = dict()
        self.errors = dict()
        #self.num_electrodes = case1.total_electrodes
        self.num_trials = len(case1.trials)

    def load_data(self):

        # range of SNR values for our plots
        self.SNR = numpy.array([snr for snr in case1.snr_values])
        self.SNR_dB = numpy.log10(self.SNR)

        # read true dipole locations
        for dp in case1.dipoles:
            filename = case1.true_dipole_dat(dp=dp)
            self.dipoles[dp] = read_dipoles(filename=case1.true_dipole_dat(dp=dp))
            assert len(self.dipoles[dp]) == 1

        # source localization errors over all combinations
        errors = dict()
        for kw in case1.every_combination():
            k = lv,dp,ec,snr,tr = case1.unpack(**kw)

            # source localization
            self.sources[k] = read_dipoles(filename=case1.sources_dat(**kw))
            assert len(self.sources[k]) == 1

            # position errors
            errors[k] = dist(self.dipoles[dp][0], self.sources[k][0])

        # now, tabulate the position errors in a form that we can more easily plot
        for lv in case1.levels:
            for dp in case1.dipoles:
                for ec in case1.electrodes:
                    position_error = numpy.zeros(shape=(self.num_trials, len(self.SNR)))
                    for i in xrange(len(self.SNR)):
                        snr = case1.snr_codes[i]
                        for j in xrange(self.num_trials):
                            tr = case1.trials[j]
                            position_error[j,i] = errors[lv,dp,ec,snr,tr]
                    self.errors[lv,dp,ec] = position_error

        return

    def plot_position_error_vs_snr(self, levels, dipoles, electrodes):

        pylab.clf()
        ax = pylab.subplot(111)
        ax.set_xscale('log')
        ax.set_yscale('log')
        pylab.xlim([0.8 * self.SNR[0], 1.25 * self.SNR[-1]])

        for lv in levels:
            for dp in dipoles:
                for ec in electrodes:
                    kw = dict(lv=lv, dp=dp, ec=ec)
                    label = '{lv} {dp} {ec}'.format(**kw)
                    position_error = 1000 * self.errors[lv,dp,ec]   # note conversion to mm
                    err = numpy.mean(position_error, axis=0)
                    err_std = numpy.std(position_error, axis=0)
                    pylab.errorbar(self.SNR, err, yerr=err_std, label=label)


        pylab.title('Spherical Head Model(radius 85 mm, brain radius 70 mm)')
        pylab.legend(loc='best')
        pylab.xlabel('SNR')
        pylab.ylabel('RMS error in current-source position (mm)')
        pylab.show()
        return


def main():

    levels = ['lv0']
    #levels = ['lv1']
    #levels = ['lv0', 'lv1']

    #electrodes = ['ec1']
    #electrodes = ['ec2']
    electrodes = ['ec1','ec2']

    #dipoles = ['dp1']
    dipoles = ['dp2']
    #dipoles = ['dp3']
    #dipoles = ['dp2','dp3']
    #dipoles = ['dp1','dp2', 'dp3']

    error_table = ErrorTable()
    error_table.load_data()
    error_table.plot_position_error_vs_snr(levels, dipoles, electrodes)

if __name__ == '__main__':
    main()

# EOF
