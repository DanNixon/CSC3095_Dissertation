# Runs the peak finder for all samples in the case study data

test_cases = {
    'ibnz_acid':'19387-19436',
    'benz_acid':'19282-19331',
    'lih':'21303-21342',
    'lid':'21143-21182',
    'sq_acid':'16929-16948',
    'bn_4k':'16648-16655',
    'bn_300k':'16656-16661',
    'c_4k':'16674-16679',
    'c_300k':'16719-16725',
    'prot_con':'14917-14928',
    'nd4':'14515-14529',
    'nh4':'14530-14539',
    'gl_zrbe':'22542-22575',
    'pc_zrbe':'22576-22608'
}

for name, runs in test_cases.iteritems():
    data = LoadVesuvio(Filename=runs,
                       SpectrumList='148',
                       Mode='SingleDifference',
                       InstrumentParFile='IP0004_10.par',
                       OutputWorkspace=name)

    FindPeaks(InputWorkspace=data,
              FWHM=25,
              Tolerance=40,
              BackgroundType='Quadratic',
              PeaksList='{0}_peaks'.format(data))
