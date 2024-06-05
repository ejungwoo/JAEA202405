#include "DetectorSetting.h"

ClassImp(DetectorSetting)

DetectorSetting* DetectorSetting::fInstance = nullptr;
DetectorSetting* DetectorSetting::GetDetectorSetting() {
    if (fInstance!=nullptr)
        return fInstance;
    return new DetectorSetting();
} 

DetectorSetting::DetectorSetting()
{
    fInstance = this;
    fMapDetectorType = new int*[fNumModules];
    fMapDetectorChannel = new int*[fNumModules];
    fMapDetectorReplaced = new bool*[fNumModules];
    fMapDetectorGroup = new int*[fNumModules];
    fMapDetectorRFMod = new int*[fNumModules];
    fMapDetectorRFMCh = new int*[fNumModules];
    for (int midx=0; midx<fNumModules; ++midx) {
        fMapDetectorType[midx] = new int[fNumChannels];
        fMapDetectorChannel[midx] = new int[fNumChannels];
        fMapDetectorReplaced[midx] = new bool[fNumChannels];
        fMapDetectorGroup[midx] = new int[fNumChannels];
        fMapDetectorRFMod[midx] = new int[fNumChannels];
        fMapDetectorRFMCh[midx] = new int[fNumChannels];
        for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
            fMapDetectorType[midx][iChannel] = kXX;
            fMapDetectorChannel[midx][iChannel] = -1;
            fMapDetectorReplaced[midx][iChannel] = false;
            fMapDetectorGroup[midx][iChannel]  = 0;
            fMapDetectorRFMod[midx][iChannel] = -1;
            fMapDetectorRFMCh[midx][iChannel] = -1;
        }
    }
    for (int iModule=0; iModule<20; ++iModule)
        fMapFEToModuleIndex[iModule] = -1;

    fMapDetectorToGlobalID = new int*[fNumDetectors];
    for (int det=0; det<fNumDetectors; ++det) {
        fMapDetectorToGlobalID[det] = new int[fMaxDCh];
        for (int dch=0; dch<fMaxDCh; ++dch) {
            fMapDetectorToGlobalID[det][dch] = -1;
        }
    }

    fMapGlobalIDToMIdx = new int[fNumCh];
    fMapGlobalIDToMCh = new int[fNumCh];
    for (int gid=0; gid<fNumCh; ++gid) {
        fMapGlobalIDToMIdx[gid] = -1;
        fMapGlobalIDToMCh[gid] = -1;
    }

    fMapGlobalIDToDet = new int[fNumCh];
    fMapGlobalIDToDCh = new int[fNumCh];
    for (int gid=0; gid<fNumCh; ++gid) {
        fMapGlobalIDToDet[gid] = -1;
        fMapGlobalIDToDCh[gid] = -1;
    }

    fDetectorName[kXX ] = "X";
    fDetectorName[kdE ] = "dEDetector";
    fDetectorName[kS1J] = "S1Junction";
    fDetectorName[kS3J] = "S3Junction";
    fDetectorName[kS3O] = "S3Ohmic";
    fDetectorName[kSC ] = "Scintillator";
    fDetectorName[kFC ] = "FaradayCup";

    TString name;
    string buffer;
    bool replaced;
    UShort_t module0, channelID0, FENumber, mch, dch, dES1Group, midx;

    std::ifstream mapFile(fMapFileName);
    if (mapFile.fail()) {
        cout << "Cannot open mapping file: " << fMapFileName << endl;
        return;
    }
    getline(mapFile, buffer); // skip header
    while (mapFile >> FENumber >> midx >> mch >> name >> dch >> dES1Group >> replaced)
    {
        int det = kXX;
        for (det=0; det<fNumDetectors; ++det)
            if (fDetectorName[det]==name)
                break;
        fMapDetectorType[midx][mch] = det;
        fMapDetectorChannel[midx][mch] = dch;
        fMapDetectorReplaced[midx][mch] = replaced;
        fMapDetectorGroup[midx][mch] = dES1Group;
        if (replaced) {
            mapFile >> module0 >> channelID0;
            fMapDetectorRFMod[midx][mch] = module0;
            fMapDetectorRFMCh[midx][mch] = channelID0;
        }
        if (fMapFEToModuleIndex[FENumber]==-1)
            fMapFEToModuleIndex[FENumber] = midx;
        auto gid = MIdxMCh_GID(midx, mch);
        fMapDetectorToGlobalID[det][dch] = gid;
        fMapGlobalIDToMIdx[gid] = midx;
        fMapGlobalIDToMCh[gid] = mch;
        fMapGlobalIDToDet[gid] = det;
        fMapGlobalIDToDCh[gid] = dch;
    }
    mapFile.close();

    Short_t strip, det;
    double rin, rout, amin, amax, amid, dist;
    std::ifstream s1File(fDetFileName);
    if (s1File.fail()) {
        cout << "!! Cannot open detter setting file: " << fDetFileName << endl;
        return;
    }
    while (s1File >> det >> dch >> strip >> rin >> rout >> dist >> amin >> amax >> amid) {
        if (det==1) {
            fMapS1ChToAngle1[dch] = amin;
            fMapS1ChToAngle2[dch] = amax;
            fMapS1ChToAngle[dch] = amid;
            fMapS1ChToStrip[dch] = strip;
        }
        if (det==3) {
            fMapS3ChToAngle1[dch] = amin;
            fMapS3ChToAngle2[dch] = amax;
            fMapS3ChToAngle[dch] = amid;
            fMapS3ChToStrip[dch] = strip;
        }
    }
}
