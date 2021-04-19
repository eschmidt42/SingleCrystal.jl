"""
Todos:
- last steps from `sites` and `kinds` to positions in real space and elements
"""

from ase.spacegroup import crystal
from ase import spacegroup
import numpy as np
import json
import warnings
import tqdm
import os
import itertools
from loguru import logger

class MiniSpacegroup:
    _no = None
    _symbol = None
    _setting = None
    _centrosymmetric = None
    _scaled_primitive_cell = None
    _reciprocal_cell = None
    _nsubtrans = None
    _subtrans = None
    _nsymop = None
    _rotations = None
    _translations = None

    def to_dict(self):
        return dict(
            no = self._no,
            symbol = self._symbol,
            centrosymmetric = self._centrosymmetric,
            scaled_primitive_cell = self._scaled_primitive_cell.tolist(),
            reciprocal_cell = self._reciprocal_cell.tolist(),
            nsubtrans = self._nsubtrans,
            subtrans = self._subtrans.tolist(),
            nsymop = self._nsymop,
            rotations = self._rotations.tolist(),
            translations = self._translations.tolist()
        )

def format_symbol(symbol):
    """Returns well formatted Hermann-Mauguin symbol as extected by
    the database, by correcting the case and adding missing or
    removing dublicated spaces."""
    fixed = []
    s = symbol.strip()
    s = s[0].upper() + s[1:].lower()
    for c in s:
        if c.isalpha():
            if len(fixed) and fixed[-1] == '/':
                fixed.append(c)
            else:
                fixed.append(' ' + c + ' ')
        elif c.isspace():
            fixed.append(' ')
        elif c.isdigit():
            fixed.append(c)
        elif c == '-':
            fixed.append(' ' + c)
        elif c == '/':
            fixed.append(c)
    s = ''.join(fixed).strip()
    return ' '.join(s.split())

class SpacegroupError(Exception):
    """Base exception for the spacegroup module."""
    pass


class SpacegroupNotFoundError(SpacegroupError):
    """Raised when given space group cannot be found in data base."""
    pass


class SpacegroupValueError(SpacegroupError):
    """Raised when arguments have invalid value."""
    pass

def _skip_to_blank(f, spacegroup, setting):
    """Read lines from f until a blank line is encountered."""
    while True:
        line = f.readline()
        if not line:
            raise SpacegroupNotFoundError(
                'invalid spacegroup `%s`, setting `%s` not found in data base'
                % (spacegroup, setting))
        if not line.strip():
            break


def _skip_to_nonblank(f, spacegroup, setting):
    """Read lines from f until a nonblank line not starting with a
    hash (#) is encountered and returns this and the next line."""
    while True:
        line1 = f.readline()
        if not line1:
            raise SpacegroupNotFoundError(
                'invalid spacegroup %s, setting %i not found in data base' %
                (spacegroup, setting))
        line1.strip()
        if line1 and not line1.startswith('#'):
            line2 = f.readline()
            break
    return line1, line2

def _read_datafile_entry(spg, no, symbol, setting, f):
    """Read space group data from f to spg."""

    floats = {'0.0': 0.0, '1.0': 1.0, '0': 0.0, '1': 1.0, '-1': -1.0}
    for n, d in [(1, 2), (1, 3), (2, 3), (1, 4), (3, 4), (1, 6), (5, 6)]:
        floats['{0}/{1}'.format(n, d)] = n / d
        floats['-{0}/{1}'.format(n, d)] = -n / d

    spg._no = no
    spg._symbol = symbol.strip()
    spg._setting = setting
    spg._centrosymmetric = bool(int(f.readline().split()[1]))
    # primitive vectors
    f.readline()
    spg._scaled_primitive_cell = np.array(
        [[float(floats.get(s, s)) for s in f.readline().split()]
         for i in range(3)],
        dtype=float)
    # primitive reciprocal vectors
    f.readline()
    spg._reciprocal_cell = np.array([[int(i) for i in f.readline().split()]
                                     for i in range(3)],
                                    dtype=int)
    # subtranslations
    spg._nsubtrans = int(f.readline().split()[0])
    spg._subtrans = np.array(
        [[float(floats.get(t, t)) for t in f.readline().split()]
         for i in range(spg._nsubtrans)],
        dtype=float)
    # symmetry operations
    nsym = int(f.readline().split()[0])
    symop = np.array([[float(floats.get(s, s)) for s in f.readline().split()]
                      for i in range(nsym)],
                     dtype=float)
    spg._nsymop = nsym
    spg._rotations = np.array(symop[:, :9].reshape((nsym, 3, 3)), dtype=int)
    spg._translations = symop[:, 9:]

def _read_datafile(spg, spacegroup, setting, f):
    if isinstance(spacegroup, int):
        pass
    elif isinstance(spacegroup, str):
        spacegroup = ' '.join(spacegroup.strip().split())
        compact_spacegroup = ''.join(spacegroup.split())
    else:
        raise SpacegroupValueError('`spacegroup` must be of type int or str')
    while True:
        line1, line2 = _skip_to_nonblank(f, spacegroup, setting)
        _no, _symbol = line1.strip().split(None, 1)
        _symbol = format_symbol(_symbol)
        compact_symbol = ''.join(_symbol.split())
        _setting = int(line2.strip().split()[1])
        _no = int(_no)
        if ((isinstance(spacegroup, int) and _no == spacegroup
             and _setting == setting)
                or (isinstance(spacegroup, str)
                    and compact_symbol == compact_spacegroup) and
            (setting is None or _setting == setting)):
            _read_datafile_entry(spg, _no, _symbol, _setting, f)
            break
        else:
            _skip_to_blank(f, spacegroup, setting)

def dat2json(json_fname='../src/spacegroup.json', dat_fname='../src/spacegroup.dat'):
    if os.path.exists(json_fname):
        logger.info(f'{json_fname} already exists.')
        return
    nums = range(1,231)
    settings = [1,2]
    spgs = {}

    for (num, setting) in tqdm.tqdm(itertools.product(nums, settings),total=len(nums)*len(settings),desc='Spacegroup+setting'):
    
        spg = MiniSpacegroup()
        with open(dat_fname, 'r') as f:
            try:
                _read_datafile(spg, num, setting, f)
            except SpacegroupNotFoundError:
                continue
            except: 
                logger.exception('Unexpected exception occurred.')
        
        spgs[f'{num}: {setting}'] = spg.to_dict()

    assert len(spgs) == 274
    logger.info(f'Writing {len(spgs)} spacegroup-setting combination to {json_fname}')
    with open(json_fname, 'w') as f:
        json.dump(spgs, f)
    
def get_symop(spg):
    """Returns all symmetry operations (including inversions and
    subtranslations) as a sequence of (rotation, translation)
    tuples."""
    symop = []
    parities = [1]
    if spg['centrosymmetric']:
        parities.append(-1)
    for parity in parities:
        for subtrans in spg['subtrans']:
            for rot, trans in zip(spg['rotations'], spg['translations']):
                newtrans = np.mod(np.array(trans) + np.array(subtrans), 1)
                symop.append((parity * np.array(rot), newtrans))
    return symop

def equivalent_sites(spg, scaled_positions,
                     onduplicates='error',
                     symprec=1e-3):
    kinds = []
    sites = []

    scaled = np.array(scaled_positions, ndmin=2)
    symops = get_symop(spg)

    for kind, pos in enumerate(scaled):
        for rot, trans in symops:
            site = np.mod(np.dot(rot, pos) + trans, 1.)
            if not sites:
                sites.append(site)
                kinds.append(kind)
                continue
            t = site - sites
            mask = np.all(
                (abs(t) < symprec) | (abs(abs(t) - 1.0) < symprec), axis=1)
            if np.any(mask):
                inds = np.argwhere(mask).flatten()
                for ind in inds:
                    # then we would just add the same thing again -> skip
                    if kinds[ind] == kind:
                        pass
                    elif onduplicates == 'keep':
                        pass
                    elif onduplicates == 'replace':
                        kinds[ind] = kind
                    elif onduplicates == 'warn':
                        warnings.warn('scaled_positions %d and %d '
                                        'are equivalent' %
                                        (kinds[ind], kind))
                    elif onduplicates == 'error':
                        raise SpacegroupValueError(
                            'scaled_positions %d and %d are equivalent' %
                            (kinds[ind], kind))
                    else:
                        raise SpacegroupValueError(
                            'Argument "onduplicates" must be one of: '
                            '"keep", "replace", "warn" or "error".')
            else:
                sites.append(site)
                kinds.append(kind)

    return np.array(sites), kinds

if __name__ == '__main__':
    json_fname='spacegroup.json'
    #dat2json(json_fname=json_fname)
    
    spgs = json.load(open(json_fname,'r'))
    
    num = 229
    setting = 1
    spg = spgs[f'{num}: {setting}']
    print(f'\nspg:\n{spg}')
    scaled_positions = [(0,0,0), ]
    sites, kinds = equivalent_sites(spg, scaled_positions,
                     onduplicates='error',
                     symprec=1e-3)
    print(f'\nsites \n{sites}\n\nkinds: \n{kinds}')
    cell = np.diag([1,1,1])*2.1
    print('cell:\n',cell)
    print('positions:\n', np.dot(sites,cell))

    num = 225
    setting = 1
    spg = spgs[f'{num}: {setting}']
    print(f'\nspg:\n{spg}')
    scaled_positions = [(0,0,0), (.5,.5,.5)]
    sites, kinds = equivalent_sites(spg, scaled_positions,
                     onduplicates='error',
                     symprec=1e-3)
    print(f'\nsites \n{sites}\n\nkinds: \n{kinds}')
    cell = np.diag([1,1,1])*2.1
    print('cell:\n',cell)
    print('positions:\n', np.dot(sites,cell))