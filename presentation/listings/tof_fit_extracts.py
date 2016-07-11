mass1 = {'value':7,  'function':'Gaussian', 'width':[4,7,10]}
mass2 = {'value':9,  'function':'Gaussian', 'width':[4,7,20]}
mass3 = {'value':19, 'function':'Gaussian', 'width':[5,10,20]}
mass4 = {'value':27, 'function':'Gaussian', 'width':13.42554}
flags['masses'] = [mass1, mass2, mass3, mass4]

flags['intenity_constraints'] = {
        [1, 0, -0.1,  0],
        [0, 1, -0.66, 0]
    }

flags['background'] = {'function':'Polynomial', 'order':2}
