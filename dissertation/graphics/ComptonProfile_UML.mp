input metauml;

beginfig(1);

iClass.iAttributeStack.iPict.iFont.scale := 0.8;
iClass.iMethodStack.iPict.iFont.scale := 0.8;


EClass.IFunction(iClassNameOnly)("IFunction")()();

EClass.ParamFunction(iClassNameOnly)("ParamFunction")()();
EClass.IFunctionOneD(iClassNameOnly)("IFunction1D")()();

EClass.ComptonProfile(iClassNameOnly)("ComptonProfile")()();

EClass.GaussianComptonProfile(iClassNameOnly)("GaussianComptonProfile")()();
EClass.GramCharlierComptonProfile(iClassNameOnly)("GramCharlierComptonProfile")()();
EClass.MultivariateGaussianComptonProfile(iClassNameOnly)("MultivariateGaussianComptonProfile")()();


Group.A(
    ParamFunction,
    IFunctionOneD
);
leftToRight(25)(
    ParamFunction,
    IFunctionOneD
);

Group.B(
    GaussianComptonProfile,
    GramCharlierComptonProfile,
    MultivariateGaussianComptonProfile
);
leftToRight(25)(
    GaussianComptonProfile,
    GramCharlierComptonProfile,
    MultivariateGaussianComptonProfile
);

topToBottom(25)(
    IFunction,
    A,
    ComptonProfile,
    B
);


drawObjects(
    IFunction,
    A,
    ComptonProfile,
    B
);


link(inheritance)(ParamFunction.n -- IFunction.s);
link(inheritance)(IFunctionOneD.n -- IFunction.s);

link(inheritance)(ComptonProfile.n -- ParamFunction.s);
link(inheritance)(ComptonProfile.n -- IFunctionOneD.s);

link(inheritance)(GaussianComptonProfile.n -- ComptonProfile.s);
link(inheritance)(GramCharlierComptonProfile.n -- ComptonProfile.s);
link(inheritance)(MultivariateGaussianComptonProfile.n -- ComptonProfile.s);

endfig;

end
