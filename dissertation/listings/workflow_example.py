from vesuvio.workflow import fit_tof

## Standard flags to modify processing

runs = "15039-15045"

flags = dict()

flags['fit_mode'] = 'spectra'
flags['spectra'] = '143-150'
flags['bin_parameters'] = None

mass1 = {'value': 1.0079, 'function': 'GramCharlier', 'width': [2, 5, 7],
         'hermite_coeffs': [1,0,0], 'k_free': 0, 'sears_flag': 1}
mass2 = {'value': 16.0, 'function': 'Gaussian', 'width': 10}
mass3 = {'value': 27.0, 'function': 'Gaussian', 'width': 13}
mass4 = {'value': 133.0, 'function': 'Gaussian', 'width': 30}
flags['masses'] = [mass1, mass2, mass3, mass4]

flags['intensity_constraints'] = list([0, 1, 0, -4])

flags['background'] = {'function': 'Polynomial', 'order': 2}

## Corrections flags

flags['output_verbose_corrections'] = True

flags['container_runs'] = None
flags['fixed_container_scaling'] = None

flags['gamma_correct'] = True
flags['fixed_gamma_scaling'] = None

flags['ms_flags'] = dict()
flags['ms_flags']['SampleWidth'] = 10.0
flags['ms_flags']['SampleHeight'] = 10.0
flags['ms_flags']['SampleDepth'] = 0.5
flags['ms_flags']['SampleDensity'] = 241

# Optional parameters (default values are given)
# flags['ms_flags']['Seed'] = 123456789
# flags['ms_flags']['NumScatters'] = 3
# flags['ms_flags']['NumRuns'] = 10
# flags['ms_flags']['NumEvents'] = 50000
# flags['ms_flags']['SmoothNeighbours'] = 3
# flags['ms_flags']['BeamRadius'] = 2.5

## Advanced flags

flags['ip_file'] = 'IP0004_10.par'
flags['diff_mode'] = 'single'

flags['max_fit_iterations'] = 5000
flags['fit_minimizer'] = 'Levenberg-Marquardt,AbsError=1e-08,RelError=1e-08'

flags['iterations'] = 1
flags['convergence'] = None

## Run fit
fit_tof(runs, flags, flags['iterations'], flags['convergence'])
