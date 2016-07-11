input metauml;

beginfig(1);

EClass.raw_file(iClassNameOnly)("Raw data file (.raw)")()();
Activity.LoadRaw("LoadRaw");
EClass.MatrixWorkspace(iClassNameOnly)("MatrixWorkspace")()();
EClass.BinParams(iClassNameOnly)("BinParameters (0.5,0.1,12)")()();
Activity.Rebin("Rebin");
EClass.MatrixWorkspace_b(iClassNameOnly)("MatrixWorkspace")()();
Activity.SaveNexusProcessed("SaveNexusProcessed");
EClass.nexus_fle(iClassNameOnly)("NeXus data file (.nxs)")()();


topToBottom(15)(
    raw_file,
    LoadRaw,
    MatrixWorkspace,
    Rebin,
    MatrixWorkspace_b,
    SaveNexusProcessed,
    nexus_fle);

leftToRight(10)(
    MatrixWorkspace,
    BinParams);


drawObjects(
    raw_file,
    LoadRaw,
    MatrixWorkspace,
    BinParams,
    Rebin,
    MatrixWorkspace_b,
    SaveNexusProcessed,
    nexus_fle);


clink(transition)(raw_file, LoadRaw);
clink(transition)(LoadRaw, MatrixWorkspace);
clink(transition)(MatrixWorkspace, Rebin);
clink(transition)(BinParams, Rebin);
clink(transition)(Rebin, MatrixWorkspace_b);
clink(transition)(MatrixWorkspace_b, SaveNexusProcessed);
clink(transition)(SaveNexusProcessed, nexus_fle);

endfig;

end
