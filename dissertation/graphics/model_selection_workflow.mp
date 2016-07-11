input metauml;

beginfig(1);
  Begin.b;
  Activity.loadData("Load data");
  Activity.findPeakAllSpec("Find Peaks in all Spectra");
  Activity.findPeakInSpec("Find Peaks in each Spectrum");
  Activity.averagePeakPositions("Create models for found peaks");
  Activity.fitModelToSpectrum("Fit model to peak");
  Activity.runFABADAOnModel("Run FABADA on model");
  End.e;


  topToBottom(25)(
      b,
      loadData,
      findPeakAllSpec,
      findPeakInSpec,
      averagePeakPositions,
      fitModelToSpectrum,
      runFABADAOnModel,
      e
  );


  drawObjects(
      b,
      loadData,
      findPeakAllSpec,
      findPeakInSpec,
      averagePeakPositions,
      fitModelToSpectrum,
      runFABADAOnModel,
      e
  );


  clink(transition)(b, loadData);
  clink(transition)(loadData, findPeakAllSpec);
  clink(transition)(findPeakAllSpec, findPeakInSpec);
  clink(transition)(findPeakInSpec, averagePeakPositions);
  clink(transition)(averagePeakPositions, fitModelToSpectrum);
  clink(transition)(fitModelToSpectrum, runFABADAOnModel);
  clink(transition)(runFABADAOnModel, e);

endfig;

end
