import mantid.simpleapi as ms
from mantid.api import *
from mantid.kernel import *
from vesuvio.commands import load_and_crop_data, _create_profile_strs_and_mass_list, _create_background_str
from vesuvio.base import VesuvioBase, TableWorkspaceDictionaryFacade
import math
import itertools
import sys
import numpy as np
import scipy

#==============================================================================

# Generated from: http://www.lfd.uci.edu/~gohlke/code/elements.py
MASSES = np.array([ 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107,
    14.0067, 15.9994, 18.9984032, 20.1797, 22.98977, 24.305, 26.981538,
    28.0855, 30.973761, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.95591,
    47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332, 58.6934, 63.546,
    65.409, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.798, 85.4678, 87.62,
    88.90585, 91.224, 92.90638, 95.94, 97.907216, 101.07, 102.9055, 106.42,
    107.8682, 112.411, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293,
    132.90545, 137.327, 138.9055, 140.116, 140.90765, 144.24, 144.912744,
    150.36, 151.964, 157.25, 158.92534, 162.5, 164.93032, 167.259, 168.93421,
    173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217,
    195.078, 196.96655, 200.59, 204.3833, 207.2, 208.98038, 208.982416,
    209.9871, 222.0176, 223.0197307, 226.025403, 227.027747, 232.0381,
    231.03588, 238.02891, 237.048167, 244.064198, 243.061373, 247.070347,
    247.070299, 251.07958, 252.08297, 257.095099, 258.098425, 259.10102,
    262.10969, 261.10875, 262.11415, 266.12193, 264.12473, 269.13411,
    268.13882])

FINAL_ENERGY = 4897 # meV
FINAL_VELOCITY = math.sqrt(FINAL_ENERGY / (5.2276 * scipy.constants.micro))
NEUTRON_MASS_AMU = scipy.constants.value("neutron mass in u")

#==============================================================================

class VesuvioModelSelection(VesuvioBase):

    _instrument_params = None
    _mass_tof_positions = None

#------------------------------------------------------------------------------

    def summary(self):
        return 'Select the most probable fitting model for VESUVIO data.'

    def category(self):
        return "Inelastic\\Indirect\\Vesuvio"

#------------------------------------------------------------------------------

    def PyInit(self):
        self.declareProperty(name='Runs', defaultValue='',
                             doc='Sample run range')

        self.declareProperty(name='FitMode', defaultValue='spectra',
                             validator=StringListValidator(['spectra', 'bank']),
                             doc='Fitting mode')

        self.declareProperty(name='Spectra', defaultValue='',
                             doc='Spectra range or bank')

        self.declareProperty(name='DiffMode', defaultValue='single',
                             validator=StringListValidator(['single', 'double']),
                             doc='Differencing mode')

        self.declareProperty(name='BinParameters', defaultValue='',
                             doc='Input binning parameters')

        self.declareProperty(name='FWHM', defaultValue=25,
                             doc='Peak finder FWHM')

        self.declareProperty(FileProperty('IPFile', '', action=FileAction.OptionalLoad,
                                          extensions=['dat', 'par']),
                             doc='Instrument parameter file')

#------------------------------------------------------------------------------

    def validateInputs(self):
        issues = dict()

        for p_name in ['Runs', 'Spectra']:
            if self.getPropertyValue(p_name) == '':
                issues[p_name] = 'Property must have a single value or range'

        return issues

#------------------------------------------------------------------------------

    def PyExec(self):
        self._mass_tof_positions = {}

        # Load instrument and parameters
        ip_file = self.getPropertyValue('IPFile');

        self._load_instrument_params(ip_file)

        bin_params = self.getPropertyValue('BinParameters')
        if bin_params == "":
            bin_params = None

        sample_data = load_and_crop_data(self.getPropertyValue('Runs'),
                                         self.getPropertyValue('Spectra'),
                                         ip_file,
                                         self.getPropertyValue('DiffMode'),
                                         self.getPropertyValue('FitMode'),
                                         bin_params)

        # Initial peak finding
        summed = ms.SumSpectra(InputWorkspace=sample_data)
        initial_peaks = self._find_peaks(summed, FWHM=self.getPropertyValue('FWHM'), Tolerance=40)
        ms.DeleteWorkspace(summed)

        # Per spectra peak finding
        peaks = []
        for p in initial_peaks:
            detected = self._find_peaks(sample_data, FWHM=25, Tolerance=8, PeakPositions=[p[1]])

            avg_peak = np.zeros(3)
            i = 0
            for dp in detected:
                spec, tof, width, intensity = dp

                if spec not in self._mass_tof_positions:
                    self._cache_mass_tof(spec)

                hwhm = width * 0.5
                tof_range = (tof - hwhm, tof + hwhm)

                mass_positions = self._mass_tof_positions[spec]
                in_fwhm = np.where(np.logical_and(tof_range[0] <= mass_positions, mass_positions <= tof_range[1]))[0]

                if len(in_fwhm) == 0:
                    continue

                avg_peak[0] += MASSES[in_fwhm[0]]
                avg_peak[1] += MASSES[in_fwhm[-1]]
                avg_peak[2] += intensity
                i += 1

            if i > 0:
                avg_peak /= i
                peaks.append(avg_peak)

        possible_peaks = (0, 3)

        # Model generation
        models = []
        for model in itertools.product(range(possible_peaks[0], possible_peaks[1] + 1), repeat=len(peaks)):
            p = []
            for i in range(len(peaks)):
                p.append((peaks[i], model[i]))
            models.append(self._generate_model(p))

        # Convert TOF to seconds
        sample_data = self._execute_child_alg("ScaleX", InputWorkspace=sample_data, OutputWorkspace=sample_data,
                                              Operation='Multiply', Factor=1e-06)

        # Fitting
        best_chi2 = sys.float_info.max
        best_model = None

        print "Num models: {}".format(len(models))
        for i, m in enumerate(models):
            try:
                chi2, model_ws = self._test_model(sample_data, m, i)
            except RuntimeError:
                continue

            print i, m, chi2

            if chi2 < best_chi2:
                best_chi2 = chi2
                best_model = model_ws

        print "Best Model (Chi2: {1:f}): {0}".format(best_model, best_chi2)

#------------------------------------------------------------------------------

    def _load_instrument_params(self, par_file):
        # Load data
        param_names = ['spectrum', 'theta', 't0', 'L0', 'L1']
        file_data = np.loadtxt(par_file, skiprows=1, usecols=[0, 2, 3, 4, 5], unpack=True)

        self._instrument_params = {}
        for name, column in zip(param_names, file_data):
            self._instrument_params[name] = column

        # Convert theta to radians
        self._instrument_params['theta'] = np.radians(self._instrument_params['theta'])

#------------------------------------------------------------------------------

    def _find_peaks(self, sample_data, **kwargs):
        peak_table = ms.FindPeaks(InputWorkspace=sample_data,
                                  PeaksList='{0}_peaks'.format(sample_data),
                                  PeakFunction='Gaussian',
                                  **kwargs)

        detected_peaks = []
        for row in range(peak_table.rowCount()):
            ws_idx = peak_table.cell('spectrum', row)
            spectrum = sample_data.getSpectrum(ws_idx).getSpectrumNo()
            centre = peak_table.cell('centre', row)
            width = peak_table.cell('width', row)
            intensity = peak_table.cell('height', row)
            detected_peaks.append((spectrum, centre, width, intensity))

        return detected_peaks

#------------------------------------------------------------------------------

    def _cache_mass_tof(self, spec):
        spec_idx = spec - 1

        t_0 = self._instrument_params['t0'][spec_idx] * scipy.constants.micro
        l_0 = self._instrument_params['L0'][spec_idx]
        l_1 = self._instrument_params['L1'][spec_idx]
        theta = self._instrument_params['theta'][spec_idx]

        r_t = (np.cos(theta) + np.sqrt((MASSES / NEUTRON_MASS_AMU)**2 - np.sin(theta)**2)) / ((MASSES / NEUTRON_MASS_AMU) + 1)
        tof = (((l_0 * r_t) + l_1) / FINAL_VELOCITY) / scipy.constants.micro

        self. _mass_tof_positions[spec] = tof

#------------------------------------------------------------------------------

    def _generate_model(self, peaks):
        model = {'function_str': 'composite=CompositeFunction,NumDeriv=1;'}
        constraints = []

        i = 0
        for p in peaks:
            for n in range(p[1]):
                mass_range = p[0][:2]
                intensity = p[0][2]
                half_mass = mass_range[0] + (mass_range[1] - mass_range[0]) * 0.5
                model['function_str'] += 'name=GaussianComptonProfile,Mass={0:f},Intensity={1:f};'.format(half_mass, intensity)
                constraints.append('f{0}.Intensity > 0.0, f{0}.Mass > {1:f}, f{0}.Mass < {2:f}, f{0}.Width > 0'.format(i, mass_range[0], mass_range[1]))
                i += 1

        model['function_str'] += 'name=Polynomial,n=2;'
        model['constraints_str'] = ','.join(constraints)

        return model

#------------------------------------------------------------------------------

    def _test_model(self, sample_data, model, model_idx):
        workspaces = []
        chi2 = []

        for ws_idx in range(sample_data.getNumberHistograms()):
            name_prefix = '{0}_{1}_'.format(sample_data, model_idx)

            # Set names
            fit_ws_name = '{0}{1}_fit'.format(name_prefix, ws_idx)
            params_name = '{0}{1}_params'.format(name_prefix, ws_idx)
            pdf_name = '{0}{1}_pdf'.format(name_prefix, ws_idx)

            # Run initial fit
            outputs = self._execute_child_alg("Fit",
                                              Function=model['function_str'],
                                              InputWorkspace=sample_data,
                                              WorkspaceIndex=ws_idx,
                                              Constraints=model['constraints_str'],
                                              CreateOutput=True,
                                              OutputCompositeMembers=True,
                                              MaxIterations=100000000,
                                              Minimizer='Levenberg-Marquardt')

            reduced_chi_square, params, fitted_data = outputs[1], outputs[3], outputs[4]
            refined_function_str = self._update_function_from_params(model['function_str'], params)

            # Run FABADA fit
            minimizer_str = 'FABADA,PDF={0}'.format(pdf_name)
            outputs = self._execute_child_alg("Fit",
                                              Function=model['function_str'],
                                              InputWorkspace=sample_data,
                                              WorkspaceIndex=ws_idx,
                                              Constraints=model['constraints_str'],
                                              CreateOutput=True,
                                              OutputCompositeMembers=True,
                                              MaxIterations=100000000,
                                              Minimizer=minimizer_str)

            reduced_chi_square, pdf, params, fitted_data = outputs[1], outputs[2], outputs[5], outputs[6]
            chi2.append(reduced_chi_square)

            # Convert fitted TOF to micro seconds
            fitted_data = self._execute_child_alg("ScaleX", InputWorkspace=fitted_data,
                                                  OutputWorkspace=fitted_data,
                                                  Operation='Multiply', Factor=1e06)

            # Store results
            AnalysisDataService.add(fit_ws_name, fitted_data)
            AnalysisDataService.add(params_name, params)
            workspaces.append(fit_ws_name)
            workspaces.append(params_name)
            if 'FABADA' in minimizer_str:
                AnalysisDataService.add(pdf_name, pdf)
                workspaces.append(pdf_name)

        group_name = '{0}model'.format(name_prefix)
        ms.GroupWorkspaces(InputWorkspaces=workspaces,
                           OutputWorkspace=group_name)

        chi2 = np.average(np.array(chi2))

        return chi2, group_name

#------------------------------------------------------------------------------

    def _update_function_from_params(self, function_str, params):
        params_dict = TableWorkspaceDictionaryFacade(params)

        functions = function_str.split(';')
        new_functions = [functions[0]]

        idx = 0;
        for func in functions[1:]:
            new_params = []

            for param in func.split(','):
                if param.split('=')[0] == 'name':
                    new_params.append(param)

            param_prefix = 'f{0}.'.format(idx)
            for name, value in params_dict.items():
                if param_prefix in name:
                    new_params.append('{0}={1}'.format(name, value))

            new_functions.append(','.join(new_params))
            idx += 1

        new_function_str = ';'.join(new_functions)

        return new_function_str

#==============================================================================

# Register algorithm with Mantid
AlgorithmFactory.subscribe(VesuvioModelSelection)
