input metauml;

beginfig(1);
  Begin.b;
  Activity.loadData("Load data");
  Activity.initialFitting("Initial fit");
  Fork.correctionsFork("h", 50);
  Activity.gammaBackground("Gamma Background");
  Activity.multipleScattering("Multiple Scattering");
  Activity.canSubtract("Container Subtraction");
  Fork.correctionsJoin("h", 50);
  Activity.fitCorrections("Fit Correction Scale Factors");
  Activity.subtractCorrections("Subtract Scaled Corrections");
  Activity.finalFitting("Final Fitting");
  End.e;


  leftToRight(20)(
      gammaBackground,
      multipleScattering,
      canSubtract
  );
  Group.corrections(
      gammaBackground,
      multipleScattering,
      canSubtract
  );

  topToBottom(25)(
      b,
      loadData,
      initialFitting,
      correctionsFork,
      corrections,
      correctionsJoin,
      fitCorrections,
      subtractCorrections,
      finalFitting,
      e
  );


  drawObjects(
      b,
      loadData,
      initialFitting,
      correctionsFork,
      corrections,
      correctionsJoin,
      fitCorrections,
      subtractCorrections,
      finalFitting,
      e
  );


  clink(transition)(b, loadData);
  clink(transition)(loadData, initialFitting);
  clink(transition)(initialFitting, correctionsFork);
  clink(transition)(correctionsFork, gammaBackground);
  clink(transition)(correctionsFork, multipleScattering);
  clink(transition)(correctionsFork, canSubtract);
  clink(transition)(gammaBackground, correctionsJoin);
  clink(transition)(multipleScattering, correctionsJoin);
  clink(transition)(canSubtract, correctionsJoin);
  clink(transition)(correctionsJoin, fitCorrections);
  clink(transition)(fitCorrections, subtractCorrections);
  clink(transition)(subtractCorrections, finalFitting);
  clink(transition)(finalFitting, e);

endfig;

end
