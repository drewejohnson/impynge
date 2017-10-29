"""

"""

import logging
import argparse

import numpy
from matplotlib import pyplot


# Load default settings for homework
__choices__ = {
    'graphite': {
        'path': 'c0.txt',
        'mass': 12,
        'density': 2.23,
        'cutoff': 12000,
        'disp': 25
    },
    'iron': {
        'path': 'fe56.txt',
        'density': 7.874,
        'mass': 56.0,
        'cutoff': 56000.0,
        'disp': 25
    }
}

MASS_NEUTRON = 1.008


def computeMassLambda(m: (float, int), M: (float, int)):
    """
    Compute the mass fraction 4*m*M/(m+M)^2."""
    return 4 * m * M / (m + M) ** 2.0


def retrieveXS(filePath, evMin=None, evMax=None):
    """Open an ENDF file and return the scattering XS"""
    logging.info('Retrieving scattering cross sections from file {}'
                .format(filePath))
    energies = []
    crossSections = []
    with open(filePath) as fp:
        line = fp.readline()
        while line[0] == '#':
            line = fp.readline()
        while line != '' and '#END' not in line:
            ev, xs = [float(xx) for xx in line.split()[:2]]
            energies.append(ev)
            crossSections.append(xs)
            line = fp.readline()
    logging.info('Done')
    energies = numpy.array(energies)
    crossSections = numpy.array(crossSections)
    bounds = energies.min(), energies.max()
    if evMin is None:
        evMin = bounds[0]
    else:
        if bounds[0] > evMin:
            logging.warning('Could not find requested minimum energy '
                           '{:.4E} eV in cross section file {}. '
                           'Using minimum found: {:.4E} eV'
                           .format(evMin, filePath, bounds[0]))
            evMin = bounds[0]

    indices = numpy.where(energies >= evMin)
    energies = energies[indices]
    crossSections = crossSections[indices]
    if evMax is None:
        evMax = bounds[1]
    else:
        if bounds[1] < evMax:
            logging.warning('Could not find requested maximum energy '
                           '{:.4E} eV in cross section file {}. '
                           'Using maximum found: {:.4E} eV'
                           .format(evMax, filePath, bounds[1]))
            evMax = bounds[1]
    indices = numpy.where(energies <= evMax)
    energies = energies[indices]
    crossSections = crossSections[indices]
    return energies, crossSections


def __parse__():
    parser = argparse.ArgumentParser()
    parser.add_argument('target', help='Name of the target',
                        choices=__choices__.keys(), type=str)
    parser.add_argument('-e', '--e-min', type=float,
                        help='Minimum neutron energy [eV]')
    parser.add_argument('-E', '--e-max', type=float,
                        help='Maximum neutron energy [eV]')
    parser.add_argument('-S', '--save', action='store_true',
                        help='Save figures to this directory')
    parser.add_argument('-N', '--no-show', action='store_true',
                        help='Do not show figures.')
    args = parser.parse_args()

    # if path is not in the arguments, assume this is the run mode
    # load from default
    defOpts = __choices__[args.target]
    args.path = defOpts['path']
    args.density = defOpts['density']
    args.mass = defOpts['mass']
    args.disp = defOpts['disp']
    args.cutoff = defOpts['cutoff']

    return args

def computeDisplacementXS(energies, scatterXS, massLambda,
                          displacementEnergy, cutoffEnergy):
    """
    Create the displacement cross-section using simple kinchen-pease.

    Parameters
    ----------
    energies: numpy.array
        energies [eV] to evaluate the cross sections at
    scatterXS: numpy.array
        scatter cross sections
    massLambda: float
        :math:`4*mM/(m+M)^2`
    displacementEnergy: float
        energy required to displace an atom from its lattice site [eV]
    cutoffEnergy: float
        Kinchin-Pease electron cutoff energy [eV]

    Returns
    -------
    dispXS: numpy.array
        Displacement XS

    """
    dispPerPKA = numpy.zeros_like(scatterXS)  # displacements per PKA
    targetEnergies = numpy.multiply(energies, massLambda) / 2.0
    denominator = 2 * displacementEnergy
    # TODO streamline this iteration
    for indx, targetKin in enumerate(targetEnergies):
        if targetKin < displacementEnergy:
            continue
        elif targetKin < 2 * displacementEnergy:
            dispPerPKA[indx] = 1.0
        elif targetKin < cutoffEnergy:
            dispPerPKA[indx] = targetKin / denominator
        elif targetKin > cutoffEnergy:
            break

    dispPerPKA[indx:] = cutoffEnergy / denominator

    return dispPerPKA * scatterXS


def __plotter__(x, y, xlabel='Incident neutron energy [eV]',
                ylabel='Cross section [b]', title=None):
    fig, ax = pyplot.subplots(1, 1)
    ax.plot(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig, ax


if __name__ == '__main__':
    args = __parse__()
    energies, scatterXS = retrieveXS(args.path, args.e_min, args.e_max)
    massLambda = computeMassLambda(MASS_NEUTRON, args.mass)
    dispXS = computeDisplacementXS(energies, scatterXS, massLambda,
                                   args.disp, args.cutoff)
    dispFig, dispAx = __plotter__(energies, dispXS)

    if not args.no_show:
        pyplot.show()
    if args.save:
        dispFig.savefig('{}_dispXS'.format(args.target))
